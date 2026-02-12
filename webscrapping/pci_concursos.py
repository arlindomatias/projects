import requests
from bs4 import BeautifulSoup

url = "https://www.pciconcursos.com.br/"
r = requests.get(url)
r.raise_for_status()

soup = BeautifulSoup(r.text, "html.parser")

links = soup.select("h3 a")

for i, link in enumerate(links[:20], 1):
    titulo = link.text.strip()
    url_noticia = link["href"]
    print(f"{i}. {titulo}")
    print(f"   {url_noticia}")