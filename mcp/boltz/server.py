from typing import Any
import httpx
from mcp.server.fastmcp import FastMCP
from starlette.applications import Starlette
from mcp.server.sse import SseServerTransport
from starlette.requests import Request
from starlette.routing import Mount, Route
from mcp.server import Server
from mcp.server.fastmcp.prompts import base
from starlette.responses import JSONResponse

from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IsoelectricPointCalculator
from Bio.Seq import Seq

import uvicorn
import re

# Initialize FastMCP server for FASTA tools (SSE)
mcp = FastMCP("FASTA")

# Constants
USER_AGENT = "FASTA-app/1.0"

VIRUS_UNIPROT_REST_API_BASE = "https://rest.uniprot.org/uniprotkb"
RCSB_DB_ENDPOINT = "https://data.rcsb.org/rest/v1/core/entry"


############################ UNIPROT Database ############################################

async def make_fasta_request(url: str) -> dict[str, Any] | None:
    """Make a request to the uniprot API with proper error handling."""
    headers = {
        "User-Agent": USER_AGENT,
        "Accept": "text/plain"
    }
    async with httpx.AsyncClient() as client:
        try:
            response = await client.get(url, headers=headers, timeout=30.0)
            response.raise_for_status()
            return response.text
        except Exception:
            return None

@mcp.tool(
        name="get_virus_protein_FASTA_format_sequence",
        description="Retrieves the amino acid sequence in FASTA format for a given viral protein using its UniProt accession code."
)
async def get_fasta_protein(uniprot_code: str) -> str:
    """
    Given a UniProt accession code for a virus protein, return its sequence in FASTA format.

    This tool retrieves the amino acid sequence of the specified viral protein using its UniProt code.
    
    Example:
        Input: "P0DTC2"
        Output: A string in FASTA format representing the spike protein of SARS-CoV-2.

    Args:
        uniprot_code: The UniProt accession code for the virus protein (e.g., "P0DTC2").

    Returns:
        A string containing the protein sequence in FASTA format.
    """
    url = f"{VIRUS_UNIPROT_REST_API_BASE}/{uniprot_code}.fasta"
    data = await make_fasta_request(url)
    return data

@mcp.tool(
    name="get_virus_protein_details",
    description="Retrieve virus protein metadata (organism, species, lineage, function) from UniProt given an accession code like 'P0DTC2'."
)
async def get_virus_protein_details(uniprot_code: str) -> dict:
    """
    Returns structured metadata about a viral protein from UniProt.
    
    Args:
        uniprot_code (str): UniProt accession code (e.g., P0DTC2).
    
    Returns:
        dict: Contains organism, scientific name, lineage, function, reference URL, RCSB structural details url etc.
    """
    url = f"{VIRUS_UNIPROT_REST_API_BASE}/{uniprot_code}.json"
    
    async with httpx.AsyncClient() as client:
        response = await client.get(url)

    if response.status_code != 200:
        return {
            "error": f"Could not retrieve data for UniProt code '{uniprot_code}'.",
            "status_code": response.status_code
        }

    entry = response.json()

    # Extract useful details
    result = {
        "accession": entry.get("primaryAccession","UNKNOWN"),
        "organism": entry.get("organism", {}).get("scientificName"),
        "lineage": " → ".join(entry.get("organism", {}).get("lineage", [])),
        "taxonomy_id": entry.get("organism", {}).get("taxonId","NONE"),
    }

    # Affected hosts
    hosts = []
    for host in entry.get("organismHosts", []):
        if "scientificName" in host:
            hosts.append(host["scientificName"]+"("+host["commonName"]+")")
    result["hosts"] = ",".join(hosts)

    # Functionality of this virus in Plaintext
    comments = entry.get("comments", [])
    for comment in comments:
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts and isinstance(texts, list):
                result["function"] = texts[0].get("value")
            break

    # Parse protein name safely
    protein_name = "Unknown"
    try:
        protein_section = entry.get("proteinDescription", {})
        recommended = protein_section.get("recommendedName") or {}
        if isinstance(recommended, list):
            recommended = recommended[0]
        full_name = recommended.get("fullName")
        if isinstance(full_name, dict):
            protein_name = full_name.get("value")
        elif isinstance(full_name, str):
            protein_name = full_name
    except Exception:
        protein_name = None

    result["virus_protein_name"] = protein_name

    # Extract sequence
    sequence = entry.get("sequence", {}).get("value")

    # Compile result
    result["virus_protein_sequence"] = sequence

    # Reference url
    reference = VIRUS_UNIPROT_REST_API_BASE+"/"+uniprot_code
    result["Reference"] = reference
    return result

@mcp.tool(
    name="analyze_sequence_properties",
    description="Analyze length, molecular weight (kDa), isoelectric point (pI), and composition of a protein sequence. Use 'get_fasta_protein' to retrieve the sequence for a UniProt ID."
)
async def analyze_sequence_properties(uniprot_code: str) -> dict:
    """
    Analyze properties of a protein sequence (raw or FASTA) for a viral protein by its UniProt accession code.

    Args:
        uniprot_code: The UniProt accession code for the virus protein (e.g., "P0DTC2").

    Returns:
        dict: {
            length: int,
            molecular_weight_kda: float,
            isoelectric_point: float,
            composition: dict
        }
    """
    data = await get_fasta_protein(uniprot_code)

    lines = data.strip().splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    clean_seq = "".join(lines).upper()

    # Validate amino acids
    if not all(res in "ACDEFGHIKLMNPQRSTVWY" for res in clean_seq):
        return {"error": "Invalid sequence. Only canonical amino acids are supported."}
    seq_obj = Seq(clean_seq)
    pI_calc = IsoelectricPointCalculator(str(seq_obj))

    return {
        "length": len(seq_obj),
        "molecular_weight_kda": round(molecular_weight(seq_obj, seq_type='protein') / 1000, 2),
        "isoelectric_point": round(pI_calc.pi(), 2),
        "composition": {aa: clean_seq.count(aa) for aa in sorted(set(clean_seq))}
    }


@mcp.tool(
    name="compare_protein_variant",
    description=(
        "Compares a mutated protein (e.g., D614G) against the reference from UniProt. "
        "Returns changes in sequence, molecular weight, pI, and composition."
    )
)
async def compare_protein_variant(uniprot_id: str, mutation: str) -> dict:
    """
    Applies a mutation like D614G to a protein sequence from UniProt and compares
    basic properties between the wildtype and the variant.

    Args:
        uniprot_id (str): UniProt accession (e.g., "P0DTC2").
        mutation (str): Mutation string in format D614G.

    Returns:
        dict: Differences in molecular weight, charge, and other properties.
    """

    try:
        fasta = await get_fasta_protein(uniprot_id)
        lines = fasta.strip().splitlines()
        if lines[0].startswith(">"):
            lines = lines[1:]
        wild_seq = "".join(lines).upper()

        match = re.match(r"([A-Z])(\d+)([A-Z])", mutation.strip().upper())
        if not match:
            return {"error": "Invalid mutation format. Use e.g., D614G."}
        orig, pos, new = match.groups()
        pos = int(pos) - 1

        if wild_seq[pos] != orig:
            return {"error": f"Reference mismatch: expected {orig} at position {pos+1}, found {wild_seq[pos]}"}

        mutated_seq = wild_seq[:pos] + new + wild_seq[pos+1:]

        # Analyze both sequences
        wild_props = analyze_sequence_properties(wild_seq)
        variant_props = analyze_sequence_properties(mutated_seq)

        return {
            "mutation": mutation,
            "wildtype": wild_props,
            "variant": variant_props,
            "position": pos + 1,
            "amino_acid_change": f"{orig} → {new}"
        }

    except Exception as e:
        return {"error": str(e)}



############################ RCSB Protein Data Bank ############################################

@mcp.tool(
    name="get_top_pdb_ids_for_uniprot",
    description="Returns up to 10 representative PDB IDs for a given UniProt protein. Useful for fetching 3D structures without flooding the client."
)
async def get_top_pdb_ids_for_uniprot(uniprot_id: str) -> list[str]:
    """
    Fetches top 10 representative PDB entries from UniProt cross-references.
    Returned PDB entries help in identifying structural details of the protein from the RCSB Database.

    Args:
        uniprot_id (str): A valid UniProt accession (e.g., 'P0DTC2').

    Returns:
        list[str]: Up to 10 unique PDB IDs linked to the protein.
    """
    url = f"{VIRUS_UNIPROT_REST_API_BASE}/{uniprot_id}.json"
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
    except Exception:
        return []

    pdb_ids = [
        xref["id"]
        for xref in data.get("uniProtKBCrossReferences", [])
        if xref.get("database") == "PDB"
    ]

    return sorted(set(pdb_ids))[:10]



@mcp.tool(
    name="get_experimental_structure_details",
    description="Fetches experimental structure metadata from RCSB PDB using a valid PDB ID (e.g., '4HHB'). Useful for grounding structure-related queries with resolution, method, and official description."
)
async def get_experimental_structure_details(pdb_id: str) -> dict:
    """
    Given a PDB ID, returns curated structure metadata from RCSB PDB.

    Args:
        pdb_id (str): The 4-character PDB ID (e.g., '4HHB').

    Returns:
        dict: Metadata including structure title, method, resolution, and download link.
    """
    url = f"{RCSB_DB_ENDPOINT}/{pdb_id}"
    
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(url, timeout=10)
            response.raise_for_status()
            entry = response.json()
    except httpx.HTTPStatusError:
        return {"error": f"No entry found for PDB ID: {pdb_id}"}
    except Exception as e:
        return {"error": f"Request failed: {str(e)}"}

    return {
        "pdb_id": pdb_id.upper(),
        "title": entry.get("struct", {}).get("title"),
        "method": entry.get("exptl", [{}])[0].get("method"),
        "resolution": next(iter(entry.get("rcsb_entry_info", {}).get("resolution_combined", [])), None),
        "deposited_atoms": entry.get("rcsb_entry_info", {}).get("deposited_atom_count"),
        "release_date": entry.get("rcsb_accession_info", {}).get("initial_release_date"),
        "keywords": entry.get("struct_keywords", {}).get("pdbx_keywords"),
        "structure_url": f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    }


@mcp.tool(
    name="get_ligand_smiles_from_uniprot",
    description=(
        "Fetches up to 10 ligands (non-polymer entities) co-crystallized with PDB structures related to a UniProt ID. "
        "Returns each ligand's SMILES, formula, and name. Useful for grounding small molecule binding partners of a protein."
    )
)
async def get_ligand_smiles_from_uniprot(uniprot_id: str) -> list[dict]:
    """
    Given a UniProt ID, returns ligand details (SMILES, formula) from top related PDB structures.

    Args:
        uniprot_id (str): A valid UniProt accession (e.g., 'P0DTC2').

    Returns:
        list[dict]: Ligand metadata including ID, name, formula, and SMILES.
    """

    pdb_ids = await(get_top_pdb_ids_for_uniprot(uniprot_id))

    ligands = []

    try:
        async with httpx.AsyncClient() as client:
            # For each PDB ID, extract ligand info from RCSB
            for pdb_id in pdb_ids:
                entry_url = f"{RCSB_DB_ENDPOINT}/{pdb_id}"

                print(entry_url)
                entry_resp = await client.get(entry_url, timeout=10)
                if entry_resp.status_code != 200:
                    continue
                entry = entry_resp.json()
                entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids", [])

                print(entity_ids)
                for eid in entity_ids:
                    ligand_url = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{eid}"
                    ligand_resp = await client.get(ligand_url, timeout=10)
                    if ligand_resp.status_code != 200:
                        continue
                    ligand_data = ligand_resp.json()

                    ligands.append(ligand_data)

        print (ligands)
        return ligands[:10]  # return top 10 ligands total

    except Exception as e:
        return [{"error": str(e)}]




# REST-style endpoint that wraps MCP tool: get_fasta_protein
async def rest_get_fasta_protein(request: Request) -> JSONResponse:
    try:
        state = request.query_params["uniprot_code"]
        result = await get_fasta_protein(state)
        return JSONResponse({"result": result})
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=400)

# REST-style endpoint that wraps MCP tool: get_virus_protein_details
async def rest_get_details_protein(request: Request) -> JSONResponse:
    try:
        state = request.query_params["uniprot_code"]
        result = await get_virus_protein_details(state)
        return JSONResponse({"result": result})
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=400)
    
# REST-style endpoint that wraps MCP tool: analyze_sequence_properties
async def rest_analyze_sequence_properties(request: Request) -> JSONResponse:
    try:
        state = request.query_params["uniprot_code"]
        result = await analyze_sequence_properties(state)
        return JSONResponse({"result": result})
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=400)
    
# REST-style endpoint that wraps MCP tool: top_pdb_ids_for_uniprot
async def rest_get_top_pdb_ids_for_uniprot(request: Request) -> JSONResponse:
    try:
        state = request.query_params["uniprot_code"]
        result = await get_top_pdb_ids_for_uniprot(state)
        return JSONResponse({"result": result})
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=400)
    

# REST-style endpoint that wraps MCP tool: top_pdb_ids_for_uniprot
async def rest_get_experimental_structure_details(request: Request) -> JSONResponse:
    try:
        state = request.query_params["pdb_id"]
        result = await get_experimental_structure_details(state)
        return JSONResponse({"result": result})
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=400)




# REST-style endpoint that wraps MCP tool: top_pdb_ids_for_uniprot
async def rest_get_ligand_smiles_from_uniprot(request: Request) -> JSONResponse:
    try:
        state = request.query_params["uniprot_code"]
        result = await get_ligand_smiles_from_uniprot(state)
        return JSONResponse({"result": result})
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=400)



@mcp.prompt()
def get_initial_prompts() -> list[base.Message]:
    return [
        base.UserMessage(
            "You are a knowledgeable bioinformatics assistant. "
            "You help users prepare molecular inputs for docking and drug discovery workflows, particularly for the Boltz system. "
            "You can fetch protein information from UniProt, download FASTA or PDB data, explain virus protein structures, and guide users on formatting inputs. "
            "When needed, you call the appropriate tools to fetch sequence data, provide metadata, or guide formatting. "
            "Always ask clarifying questions if input is ambiguous. "
            "Keep your responses concise, scientific, and user-friendly."
        )
    ]



def create_starlette_app(mcp_server: Server, *, debug: bool = False) -> Starlette:
    """Create a Starlette application that can server the provied mcp server with SSE."""
    sse = SseServerTransport("/messages/")

    async def handle_sse(request: Request) -> None:
        async with sse.connect_sse(
                request.scope,
                request.receive,
                request._send,  # noqa: SLF001
        ) as (read_stream, write_stream):
            await mcp_server.run(
                read_stream,
                write_stream,
                mcp_server.create_initialization_options(),
            )

    return Starlette(
        debug=debug,
        routes=[
            Route("/sse", endpoint=handle_sse),
            Route("/rest/get_fasta", endpoint=rest_get_fasta_protein, methods=["GET"]),
            Route("/rest/get_protein_details", endpoint=rest_get_details_protein, methods=["GET"]),
            Route("/rest/analyze_sequence_properties", endpoint=rest_analyze_sequence_properties, methods=["GET"]),
            Route("/rest/top_pdb_ids", endpoint=rest_get_top_pdb_ids_for_uniprot, methods=["GET"]),
            Route("/rest/get_experimental_structure_details", endpoint=rest_get_experimental_structure_details, methods=["GET"]),
            Route("/rest/get_ligand_smiles_from_uniprot", endpoint=rest_get_ligand_smiles_from_uniprot, methods=["GET"]),
            Mount("/messages/", app=sse.handle_post_message),
        ],
    )



if __name__ == "__main__":
    mcp_server = mcp._mcp_server  # noqa: WPS437

    import argparse
    
    parser = argparse.ArgumentParser(description='Run MCP SSE-based server')
    parser.add_argument('--host', default='0.0.0.0', help='Host to bind to')
    parser.add_argument('--port', type=int, default=8080, help='Port to listen on')
    args = parser.parse_args()

    # Bind SSE request handling to MCP server
    starlette_app = create_starlette_app(mcp_server, debug=True)

    uvicorn.run(starlette_app, host=args.host, port=args.port)