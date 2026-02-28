import requests
import csv

def fetch_disgenet_diseases():
    api_key = "73587582-a451-4aab-99e3-f3ef27a619c3"
    url = "https://www.disgenet.org/api/v1/diseases"

    headers = {
        "Accept": "application/json",
        "Authorization": f"Bearer {api_key}"
    }

    try:
        response = requests.get(url, headers=headers)

        print("Status code:", response.status_code)
        print("Response snippet:", response.text[:500])

        if response.status_code == 200:
            diseases = response.json()
            return diseases
        else:
            print(f"Failed to fetch diseases. Status code: {response.status_code}")
            print("Response:", response.text)
            return None
    except requests.exceptions.RequestException as e:
        print("Request failed:", e)
        return None
    except ValueError as ve:
        print("JSON decode failed:", ve)
        return None

if __name__ == "__main__":
    diseases = fetch_disgenet_diseases()
    if diseases:
        print(f"Fetched {len(diseases)} diseases.")
        keys = diseases[0].keys()
        with open("../data/diseases.csv", "w", newline="", encoding="utf-8") as f:
            dict_writer = csv.DictWriter(f, fieldnames=keys)
            dict_writer.writeheader()
            dict_writer.writerows(diseases)
        print("Saved diseases.csv")
    else:
        print("No data fetched.")
