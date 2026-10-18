#!/usr/bin/env python3

import requests
from datetime import datetime, timezone
from collections import Counter

# Import the API key from the .env file
from dotenv import load_dotenv
import os
load_dotenv()

# Initial settings
lat = -21.1909725
lon = -47.8261484


# Get the API key
openweather_key = os.getenv("OPENWEATHER")


# Main function to get current weather data
def now(latitude=lat, longitude=lon, api_key=openweather_key):

    # Build the URL
    base_url = "http://api.openweathermap.org/data/2.5/weather?"
    params = f"lat={latitude}&lon={longitude}&appid={api_key}&units=metric&lang=en"
    url = base_url + params

    try:
        # API request using requests (simpler than urllib)
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()

            current_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            weather = data["weather"][0]["description"]
            temperature = data["main"]["temp"]
            humidity = data["main"]["humidity"]
            feels_like = data["main"]["feels_like"]
            temp_min = data["main"]["temp_min"]
            temp_max = data["main"]["temp_max"]
            pressure = data["main"]["pressure"]
            wind_speed = data["wind"]["speed"]
            wind_dir = data["wind"]["deg"]

            # Fix: explicit timezone
            sunrise = datetime.fromtimestamp(data["sys"]["sunrise"], tz=timezone.utc).strftime("%H:%M:%S")
            sunset = datetime.fromtimestamp(data["sys"]["sunset"], tz=timezone.utc).strftime("%H:%M:%S")

            print(f"Weather: {weather}\n"
                  f"Time: {current_time}\n"
                  f"Temperature: {temperature} °C\n"
                  f"Humidity: {humidity} %\n"
                  f"Feels like: {feels_like} °C\n"
                  f"Pressure: {pressure} hPa\n"
                  f"Wind: {wind_dir}° {wind_speed} m/s\n"
                  f"Sunrise/Sunset: {sunrise} UTC / {sunset} UTC\n")
        else:
            print(f"API error: {response.status_code}, {response.reason}")
    except requests.exceptions.RequestException as e:
        print(f"Error accessing the API: {e}")
    start()


# Main function for the weather forecast
def forecast(latitude=lat, longitude=lon, api_key=openweather_key):
    # Build the URL
    base_url = "http://api.openweathermap.org/data/2.5/forecast?"
    params = f"lat={latitude}&lon={longitude}&appid={api_key}&units=metric&lang=en"
    url = base_url + params

    try:
        # API request using requests
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            forecasts = data["list"]

            processed_data = []
            for item in forecasts:
                date = item["dt_txt"][:10]  # Extract the date
                temp_min = item["main"]["temp_min"]
                temp_max = item["main"]["temp_max"]
                weather = item["weather"][0]["description"]
                processed_data.append({"date": date, "temp_min": temp_min, "temp_max": temp_max, "weather": weather})

            # Group by date
            grouped_data = {}
            for entry in processed_data:
                date = entry["date"]
                if date not in grouped_data:
                    grouped_data[date] = {
                        "temp_min": [],
                        "temp_max": [],
                        "weather": []
                    }
                grouped_data[date]["temp_min"].append(entry["temp_min"])
                grouped_data[date]["temp_max"].append(entry["temp_max"])
                grouped_data[date]["weather"].append(entry["weather"])

            summary = []
            for date, values in grouped_data.items():
                temp_min = min(values["temp_min"])
                temp_max = max(values["temp_max"])
                weather = Counter(values["weather"]).most_common(1)[0][0]
                summary.append({"date": date, "temp_min": temp_min, "temp_max": temp_max, "weather": weather})

            summary = sorted(summary, key=lambda x: x["date"])[:5]

            for day in summary:
                print(f"Date: {datetime.strptime(day['date'], '%Y-%m-%d').strftime('%d/%m/%Y')}")
                print(f"Minimum temperature: {day['temp_min']} °C")
                print(f"Maximum temperature: {day['temp_max']} °C")
                print(f"Weather: {day['weather']}")
                print("----")
        else:
            print(f"API error: {response.status_code}, {response.reason}")
    except requests.exceptions.RequestException as e:
        print(f"Error accessing the API: {e}")
    start()


def experiments():
    with open(r"C:\Users\AsRock\Google Drive\Doutorado\Resultados\Bioinformática\Scripts\experimentos.py",
              encoding="utf-8") as f:
        script = f.read()

    # Create a global scope dictionary for exec()
    global_scope = globals().copy()

    # Run the script in the current global scope
    exec(script, global_scope)


def start():
    print("\nWhat would you like to know?")
    print("1: Current weather;")
    print("2: Weather forecast;")
    print("3: Run calculations;")
    print("4: Exit")

    type_option = input("Enter the option: \n")

    if type_option == "1":
        now()
    elif type_option == "2":
        forecast()
    elif type_option == "3":
        experiments()
    elif type_option == "4":
        print("Exiting...")
    else:
        print("Option not found.")
        start()


start()