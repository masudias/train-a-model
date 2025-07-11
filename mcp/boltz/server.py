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

import uvicorn

# Initialize FastMCP server for FASTA tools (SSE)
mcp = FastMCP("FASTA")

# Constants
USER_AGENT = "FASTA-app/1.0"

VIRUS_UNIPROT_REST_API_BASE = "https://rest.uniprot.org/uniprotkb"


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

@mcp.tool(description="Retrieves the amino acid sequence in FASTA format for a given viral protein using its UniProt accession code.")
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



# REST-style endpoint that wraps MCP tool: get_fasta_protein
async def rest_get_fasta_protein(request: Request) -> JSONResponse:
    try:
        state = request.query_params["uniprot_code"]
        result = await get_fasta_protein(state)
        return JSONResponse({"result": result})
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=400)

@mcp.prompt()
def get_initial_prompts() -> list[base.Message]:
    return [
        base.SystemMessage(
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