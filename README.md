# tema1

Tema 1 APD - Paduretu David-Pavel 332CC

Implementare:
    Am avut de paralelizat functiile : rescale_image( pentru a rescala imaginile mai mari de 2048 x 2048), sample_grid si march.
    
    Am paralelizat functia rescale_image in corpul functiei date pe schelet, am modificat argumentele functiei si tipul , am adaugat long id (id ul thread ului ) si int P(nr total de thread uri ). Am facut asta pt a putea modifica for urile pt a folosi formulele de start si end din cadrul laboratorului pt a face paralelizarea. Am folosit acelasi procedeu si pentru sample_grid si march doar ca acestea au fost implementate local in corpul functiei thread_func.
    In thread_func verific initial daca trebuie facut rescale pe imagine sau nu si pe corpul fiecarui if fac algoritmii de sample_grid si march de fiecare data calculez start si end pt  fiecare for. Facem asta pt a imparti in intervale pentru fiecare thread.

    Pentru threaduri am creat o structura care contine informatiile independente de thread dar si cele specifice fiecarui thread(id) , P nr total de thread uri, toate thread urile primesc aceeasi image, contour_map, step_x, step_y, grid. Am inclus si campul scaled_image. Acolo salvam daca imaginea trebuie scalata si lucram cu scaled image in pasii de sample si march. La final verificam dimensiunile initiale ale imaginii si daca trebuie scalata, dam write in output la campul image fiindca algoritmul s a aplicat pe el, daca nu dam write la scaled_image.
