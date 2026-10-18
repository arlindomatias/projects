import webbrowser
import time
from urllib.parse import quote_plus

QUERIES = [
    "Concurso professor doutor",
    "Concurso pesquisador",
    "Concurso professor instituto federal",
    "concurso professor universidade federal",
    "concurso professor farmacologia",
    "concurso professor bioquímica",
    "concurso professor neurociência",
]

def url_noticias(query, recentes_primeiro=True, periodo=None):
    # tbm=nws -> guia Notícias
    # sbd:1   -> ordenar por data (mais recentes primeiro)
    # qdr:d/w/m -> opcional: último dia/semana/mês
    tbs = ["sbd:1"] if recentes_primeiro else []
    if periodo:
        tbs.append(f"qdr:{periodo}")
    url = f"https://www.google.com/search?q={quote_plus(query)}&tbm=nws&hl=pt-BR"
    if tbs:
        url += "&tbs=" + ",".join(tbs)
    return url

for q in QUERIES:
    webbrowser.open_new_tab(url_noticias(q))
    time.sleep(1)  # evita abrir tudo de uma vez e disparar CAPTCHA