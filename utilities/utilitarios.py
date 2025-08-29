import requests
from datetime import datetime, timezone
from collections import Counter

# Definições iniciais
lat = -21.1909725
lon = -47.8261484


# Função para obter a chave da API
def get_key(service):
    api_keys = {
        "OpenWeather": "",
    }
    return api_keys.get(service)


# Função principal para obter os dados do clima
def agora(latitude=lat, longitude=lon, api_key=None):
    if api_key is None:
        api_key = get_key("OpenWeather")

    # Construção da URL
    base_url = "http://api.openweathermap.org/data/2.5/weather?"
    params = f"lat={latitude}&lon={longitude}&appid={api_key}&units=metric&lang=pt_br"
    url = base_url + params

    try:
        # Requisição à API usando requests (mais simples que urllib)
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()

            hora = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            clima = data["weather"][0]["description"]
            temperatura = data["main"]["temp"]
            umidade = data["main"]["humidity"]
            sensacao = data["main"]["feels_like"]
            temp_min = data["main"]["temp_min"]
            temp_max = data["main"]["temp_max"]
            pressao = data["main"]["pressure"]
            vento_vel = data["wind"]["speed"]
            vento_dir = data["wind"]["deg"]

            # Correção para timezone explícito
            nascer_sol = datetime.fromtimestamp(data["sys"]["sunrise"], tz=timezone.utc).strftime("%H:%M:%S")
            por_sol = datetime.fromtimestamp(data["sys"]["sunset"], tz=timezone.utc).strftime("%H:%M:%S")

            print(f"Clima: {clima}\n"
                  f"Hora: {hora}\n"
                  f"Temperatura: {temperatura} °C\n"
                  f"Umidade: {umidade} %\n"
                  f"Sensação térmica: {sensacao} °C\n"
                  f"Pressão: {pressao} hPa\n"
                  f"Vento: {vento_dir}° {vento_vel} m/s\n"
                  f"Nascer/Pôr-do-sol: {nascer_sol} UTC / {por_sol} UTC\n")
        else:
            print(f"Erro na API: {response.status_code}, {response.reason}")
    except requests.exceptions.RequestException as e:
        print(f"Erro ao acessar a API: {e}")
    inicio()


# Função principal para previsão do tempo
def previsao(latitude=lat, longitude=lon, api_key=None):
    if api_key is None:
        api_key = get_key("OpenWeather")

    # Construção da URL
    base_url = "http://api.openweathermap.org/data/2.5/forecast?"
    params = f"lat={latitude}&lon={longitude}&appid={api_key}&units=metric&lang=pt_br"
    url = base_url + params

    try:
        # Requisição à API usando requests
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            forecasts = data["list"]

            processed_data = []
            for forecast in forecasts:
                date = forecast["dt_txt"][:10]  # Extraindo a data
                temp_min = forecast["main"]["temp_min"]
                temp_max = forecast["main"]["temp_max"]
                clima = forecast["weather"][0]["description"]
                processed_data.append({"date": date, "temp_min": temp_min, "temp_max": temp_max, "clima": clima})

            # Agrupando por data
            grouped_data = {}
            for entry in processed_data:
                date = entry["date"]
                if date not in grouped_data:
                    grouped_data[date] = {
                        "temp_min": [],
                        "temp_max": [],
                        "clima": []
                    }
                grouped_data[date]["temp_min"].append(entry["temp_min"])
                grouped_data[date]["temp_max"].append(entry["temp_max"])
                grouped_data[date]["clima"].append(entry["clima"])

            resumo = []
            for date, values in grouped_data.items():
                temp_min = min(values["temp_min"])
                temp_max = max(values["temp_max"])
                clima = Counter(values["clima"]).most_common(1)[0][0]
                resumo.append({"date": date, "temp_min": temp_min, "temp_max": temp_max, "clima": clima})

            resumo = sorted(resumo, key=lambda x: x["date"])[:5]

            for dia in resumo:
                print(f"Data: {datetime.strptime(dia['date'], '%Y-%m-%d').strftime('%d/%m/%Y')}")
                print(f"Temperatura Mínima: {dia['temp_min']} °C")
                print(f"Temperatura Máxima: {dia['temp_max']} °C")
                print(f"Clima: {dia['clima']}")
                print("----")
        else:
            print(f"Erro na API: {response.status_code}, {response.reason}")
    except requests.exceptions.RequestException as e:
        print(f"Erro ao acessar a API: {e}")
    inicio()


def experimentos():
    with open(r"C:\Users\AsRock\Google Drive\Doutorado\Resultados\Bioinformática\Scripts\experimentos.py",
              encoding="utf-8") as f:
        script = f.read()

    # Criando um dicionário de escopo global para o exec()
    global_scope = globals().copy()

    # Executa o script no escopo global atual
    exec(script, global_scope)


def inicio():
    print("\nO que deseja saber?")
    print("1: Tempo agora;")
    print("2: Previsão do tempo;")
    print("3: Fazer cálculos;")
    print("4: Sair")

    type_option = input("Digite a opção: \n")

    if type_option == "1":
        agora()
    elif type_option == "2":
        previsao()
    elif type_option == "3":
        experimentos()
    elif type_option == "4":
        print("Saindo...")
    else:
        print("Opção não encontrada.")
        inicio()


inicio()