def oblicz_srednia(liczby):
    return sum(liczby) / len(liczby)


liczby = []

print("Podaj liczby. Wpisz 'koniec', aby zakończyć.")

while True:
    wpis = input("Liczba: ")

    if wpis.lower() == "koniec":
        break

    try:
        liczby.append(float(wpis))
    except ValueError:
        print("To nie jest poprawna liczba. Spróbuj ponownie.")

if liczby:
    srednia = oblicz_srednia(liczby)
    print(f"Średnia podanych liczb wynosi: {srednia:.2f}")
else:
    print("Nie podano żadnych liczb.")