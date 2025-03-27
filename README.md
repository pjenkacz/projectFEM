Program obliczeniowy metodą elementów skończonych w C++​
--------------------------------------------------------

**Opis​**

Niniejszy program jest zaawansowanym narzędziem do przeprowadzania analiz numerycznych z wykorzystaniem metody elementów skończonych (MES). Został napisany w języku C++ i umożliwia modelowanie oraz symulację różnorodnych problemów inżynierskich poprzez dyskretyzację ciągłych domen na skończoną liczbę elementów.​

**Funkcjonalności​**

*   **Wczytywanie danych wejściowych**: Program odczytuje informacje o siatce elementów, parametrach materiałowych oraz warunkach brzegowych z plików tekstowych.​
    
*   **Obliczenia macierzy elementów**: Dla każdego elementu programu oblicza macierze sztywności (H), macierze brzegowe (Hbc), wektory obciążeń (P) oraz macierze pojemności cieplnej (C).​
    
*   **Agregacja macierzy globalnych**: Lokalne macierze elementów są integrowane w celu utworzenia globalnych macierzy układu.​
    
*   **Rozwiązywanie układu równań**: Program wykorzystuje eliminację Gaussa do rozwiązania układu równań liniowych, co pozwala na uzyskanie rozkładu temperatur w zadanym obszarze.​
    
*   **Symulacja w czasie**: Możliwość przeprowadzania analiz w dziedzinie czasu poprzez iteracyjne aktualizowanie stanu układu w kolejnych krokach czasowych.​
    

**Wykorzystane technologie i biblioteki​**

*   **Standardowa biblioteka C++**: Program korzysta z podstawowych komponentów standardowej biblioteki C++, takich jak:​
    
    *   : do obsługi wejścia i wyjścia,​
        
    *   : do dynamicznych struktur danych,​
        
    *   : do operacji na plikach,​
        
    *   : do obliczeń matematycznych,​
        
    *   : do manipulacji strumieniami,​
        
    *   : do formatowania wyjścia,​
        
    *   : do operacji na łańcuchach znaków,​
        
    *   : do algorytmów standardowych.​
        

**Struktura programu​**

*   **Pliki nagłówkowe**: Deklaracje struktur danych oraz funkcji pomocniczych znajdują się w plikach nagłówkowych, takich jak header.h.​
    
*   **Główna funkcja programu (main)**: Odpowiada za inicjalizację danych, wywoływanie odpowiednich funkcji obliczeniowych oraz zarządzanie przepływem programu.​
    
*   **Struktury danych**: Program definiuje szereg struktur, takich jak Wezly, Jakobian, Macierz, Punkt, GlobalData, Element, które reprezentują odpowiednio węzły siatki, macierze Jacobiego, macierze elementów, punkty w przestrzeni, globalne dane symulacji oraz poszczególne elementy siatki.​
    

**Konfiguracja​**

Program umożliwia konfigurację poprzez preprocesorowe definicje:​

*   #define npc 16: Określa liczbę punktów całkowania w schemacie Gaussa (np. 4, 9 lub 16).​
    
*   #define NPCBC 4: Ustala liczbę punktów całkowania na krawędziach dla warunków brzegowych (np. 2, 3 lub 4).​
    
*   #define TEST 2: Wybór zestawu danych testowych (np. 1, 2 lub 3), co pozwala na przeprowadzanie symulacji na różnych siatkach i konfiguracjach.​
    

Sposób użycia​

1.  **Przygotowanie danych wejściowych**: Utwórz plik tekstowy zawierający informacje o węzłach, elementach, parametrach materiałowych oraz warunkach brzegowych zgodnie z oczekiwanym formatem.​
    
2.  **Kompilacja programu**: Użyj kompilatora zgodnego ze standardem C++ (np. g++) do skompilowania kodu źródłowego.​
    
3.  **Uruchomienie programu**: Po skompilowaniu uruchom program, podając jako argument ścieżkę do pliku z danymi wejściowymi.​
    
4.  **Analiza wyników**: Program wygeneruje na standardowym wyjściu informacje o przebiegu symulacji, w tym minimalne i maksymalne wartości temperatur w kolejnych krokach czasowych.
    

Autor: Michał Pieniek​
