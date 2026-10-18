import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List

from dotenv import load_dotenv
from litellm import completion

load_dotenv()
gemini_key = os.getenv("GEMINI")

PASTA_CONVERSAS = Path("conversas")

SYSTEM_PROMPT = """Você é um agente de IA especializado em data science e bioinformática.
Responda de forma completa e detalhada, sem cortar explicações."""


def generate_response(messages: List[Dict]) -> str:
    """Chama o LLM e retorna a resposta."""
    response = completion(
        model="gemini/gemini-2.5-flash",
        messages=messages,
        max_tokens=10000,  # Pode ser até 65000
        api_key=gemini_key,
    )
    return response.choices[0].message.content


def salvar_conversa(messages: List[Dict]) -> Path | None:
    """Salva o histórico em .txt (ignora se só houver o system prompt)."""
    if len(messages) <= 1:
        return None

    PASTA_CONVERSAS.mkdir(exist_ok=True)
    agora = datetime.now()
    caminho = PASTA_CONVERSAS / agora.strftime("conversa_%Y%m%d_%H%M%S.txt")

    rotulos = {"system": "SISTEMA", "user": "VOCÊ", "assistant": "AGENTE"}

    with open(caminho, "w", encoding="utf-8") as f:
        f.write(f"Conversa de {agora:%d/%m/%Y %H:%M}\n")
        f.write("=" * 60 + "\n\n")
        for m in messages:
            rotulo = rotulos.get(m["role"], m["role"].upper())
            f.write(f"[{rotulo}]\n{m['content']}\n\n")
            f.write("-" * 60 + "\n\n")

    print(f"Conversa salva em: {caminho}")
    return caminho


def chat():
    messages = [{"role": "system", "content": SYSTEM_PROMPT}]
    print("Agente de Bioinformática iniciado. Digite 'sair' para encerrar.\n")

    try:
        while True:
            problema = input("Você: ").strip()

            if problema.lower() in ("sair", "exit", "quit"):
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
                messages.pop()  # remove a pergunta que falhou
    except (KeyboardInterrupt, EOFError):
        print("\nInterrompido.")
    finally:
        salvar_conversa(messages)  # salva em qualquer saída


if __name__ == "__main__":
    chat()