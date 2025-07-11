# Boltz based bio-informatics AI agent

## Overview

Agentic AI lets people control system components using plain language. It’s flexible enough that experts from different fields can connect it to their own tools—helping reduce hallucinations by ensuring the AI talks to real systems, not just guessing. On that ground, we built this project for bioinformatics researchers who want to explore drug discovery, analyze protein structures, prototype ideas quickly, and more—all through natural language commands. No need to learn every tool’s syntax or API.

For example, instead of writing code to fetch a protein or run a molecule generation tool, you can simply ask:

```bash
Generate 5 small molecules that bind to the spike protein of SARS-CoV-2.
```
Behind the scenes, the AI (i.e. LLM) connects to real tools like Boltz, Pocket2Mol, or ESMFold to get the job done—accurately and reproducibly.

Currently we are focusing only integrating LLM to explore the vast capability of the `boltz` tool. **Contributions** are welcome to help expand integration with popular tools like `Pocket2Mol`, `DeepChem`, and `RDKit` —with the goal of supporting a broader range of bioinformatics workflows.

We leverage Model Control Protocl (mcp) to develop the bio-informatics AI agent. 

## Quick Start

### Installation

Create a virtual environment using `python` version 3.12 or later. 

```bash
python -m venv mcp
```
Then activate the environment *mcp*.

```bash
source mcp/bin/activate
```

Then install the rqquired packages:

```bash
pip install -r requirements.txt
```

Lastly, run the server code:

```bash
python run boltz/server.py --host localhost --port 8080
```

### Basic Usage

Although an LLM client (e.g., `mcp-client`) is typically needed to test our MCP server, relying on an external `API_KEY` from services like OpenAI or Anthropic can slow down local development. To reduce this dependency, we've added a simple REST interface that allows us to test and develop our tools **without requiring access to a live MCP client**. This makes it easier to create and validate test cases before exposing the service to real LLMs. The REST service is simply a wrapper around the `mcp` tool. 

Use the following `curl` request to test the REST endpoint:

```bash
curl "http://localhost:8080/rest/get_fasta?uniprot_code=P09261"
```

The string *P09261* is the UniProt code for a `SARS-CoV-2` protein.
The REST service will return the protein's amino acid sequence in `FASTA` format.

The `@mcp.prompt()` and the `@mcp.tool` decorator are written in a way that allows an LLM client to clearly understand the task and invoke the appropriate tool.