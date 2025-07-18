import asyncio
import json
import os
from typing import Optional
from contextlib import AsyncExitStack

from mcp import ClientSession
from mcp.client.sse import sse_client

from openai import OpenAI
from openai.types import Completion
from dotenv import load_dotenv

load_dotenv()  # load environment variables from .env

class MCPClient:
    def __init__(self):
        # Initialize session and client objects
        self.session: Optional[ClientSession] = None
        self.exit_stack = AsyncExitStack()
        self.openai = OpenAI(
            api_key=os.getenv("OPENAI_API_KEY")
        )
        self.available_tools = []
        self.messages = []

    async def connect_to_sse_server(self, server_url: str):
        """Connect to an MCP server running with SSE transport"""
        print("Connecting to MCP SSE server...")
        # Store the context managers so they stay alive
        self._streams_context = sse_client(url=server_url)
        streams = await self._streams_context.__aenter__()
        print("Streams:", streams)  

        self._session_context = ClientSession(*streams)
        self.session: ClientSession = await self._session_context.__aenter__()

        # Initialize
        print("Initializing SSE client...")
        await self.session.initialize()
        print("Initialized SSE client")
        
        # List available tools to verify connection
        await self.get_available_tools();
        await self.get_initial_prompts();
    
    async def cleanup(self):
        """Properly clean up the session and streams"""
        if self._session_context:
            await self._session_context.__aexit__(None, None, None)
        if self._streams_context:
            await self._streams_context.__aexit__(None, None, None)

    async def get_initial_prompts(self):
        prompt = await self.session.get_prompt("get_initial_prompts")
        # Format messages for OpenAI
        messages = []
        for message in prompt.messages:
            messages.append({
                "role": message.role,
                "content": message.content.text
            })
        self.messages = messages

    async def get_available_tools(self):
        """Get available tools from the server"""
        print("Fetching available server tools...")
        response = await self.session.list_tools()
        print("Connected to MCP server with tools:", [tool.name for tool in response.tools])

        # Format tools for OpenAI
        available_tools = [
            {
                "type": 'function',
                "function": {
                    "name": tool.name,
                    "description": tool.description,
                    "parameters": tool.inputSchema,
                },
                "strict": True,
            }
            for tool in response.tools
        ]
        self.available_tools = available_tools;


    async def call_openai(self) -> str:
        """Call OpenAI with the current messages and available tools"""
        response = self.openai.chat.completions.create(
            model="gpt-4.1-nano", # Fix these configurations to tune your response. FYI the gpt-4.1-nano cheapest one.Check the available options :https://platform.openai.com/docs/pricing
            max_tokens=1000,
            messages=self.messages,
            tools=self.available_tools
        )
        return response

    async def process_openai_response(self, response: Completion) -> str:
        """Process the response from OpenAI"""
        for choice in response.choices:
            if choice.finish_reason == "tool_calls":
                # We need to include the original message and assistant response
                self.messages.append(choice.message)
                
                for tool_call in choice.message.tool_calls:
                    tool_name = tool_call.function.name
                    tool_args = json.loads(tool_call.function.arguments)

                    #print(f"\n[Calling tool {tool_name} with args {tool_args}]...")
                    result = await self.session.call_tool(tool_name, tool_args)
                    #print(f"\nTool response: waiting..")
                    self.messages.append(
                        {
                            "role": "tool",
                            "tool_call_id": tool_call.id,
                            "content": result.content,
                        }
                    )
                
                response = await self.call_openai()
                return await self.process_openai_response(response)

            elif choice.finish_reason == "stop":
                print("\nAssistant: " + choice.message.content)
                return choice.message.content

    async def process_query(self, query: str) -> str:
        """Process a query using Claude and available tools"""
        self.messages.append({
            "role": "user",
            "content": query
        })

        response = await self.call_openai()
        return await self.process_openai_response(response)

    async def chat_loop(self):
        """Run an interactive chat loop"""
        print("\nMCP Client Started!")
        print("Type your queries or 'quit' to exit.")
        
        while True:
            try:
                print("\n" + "-" * 100)
                query = input("\nQuery: ").strip()
                
                if query.lower() == 'quit':
                    break
                    
                if query:
                    await self.process_query(query)
                    
            except Exception as e:
                print(f"\nError: {str(e)}")


async def main():
    client = MCPClient()
    try:
        await client.connect_to_sse_server(server_url=os.getenv("MCP_SSE_URL"))
        await client.chat_loop()
    finally:
        await client.cleanup()


if __name__ == "__main__":
    asyncio.run(main()) 