from litellm import completion
from typing import List, Dict
import json
from datetime import datetime
from keys import Gemini

gemini_key = Gemini

SYSTEM_PROMPT = """Você é um agente de IA especializado em data science e bioinformática.
Responda de forma completa e detalhada, sem cortar explicações."""

def generate_response(messages: List[Dict]) -> str:
    """Call LLM to get response"""
    response = completion(
        model="gemini/gemini-2.5-flash",
        messages=messages,
        max_tokens=10000, # Pode ser até 65000
        api_key=gemini_key
    )
    return response.choices[0].message.content

def chat():
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT}
    ]

    print("Agente de Bioinformática iniciado. Digite 'sair' para encerrar.\n")

    while True:
        problema = input("Você: ").strip()

        if problema.lower() in ("sair", "exit", "quit"):
            salvar_conversa(messages)  # salva ao sair
            print("Encerrando agente.")
            break

        if not problema:
            continue

        messages.append({"role": "user", "content": problema})

        try:
            resposta = generate_response(messages)
            print(f"\nAgente: {resposta}\n")
            messages.append({"role": "assistant", "content": resposta})

        except Exception as e:
            print(f"Erro ao chamar a API: {e}")


if __name__ == "__main__":
    chat()