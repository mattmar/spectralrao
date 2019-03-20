#spectralIndexes(...)
Work in progress

By Martina:

Il codice è costituito da una funzione principale e 6 subordinate, utilizzate per realizzare il calcolo (3 per nc.cores = 1 [S], 3 per nc.cores > 1 [P]). Sono predisposti 3 "mode":
-"single": un solo valore di alpha come input e una sola matrice con l'indice corrispondente in output;
-"iterative": alpha è un vettore con valore di inizio e di fine e come output una lista di matrici con gli indici corrispondenti a tutti i valori interi di alpha compresi gli estremi;
-"sequential": alpha è un vettore di lunghezza almeno 2 e come output una lista di matrici con indici corrispondenti ai valori di alpha assegnati.
Sono presenti 2 "type": "Rényi" per calcolare l'indice di Rényi e "Hill" per quello di Hill.
È presente un booleano "integer" per imporre solamente valori di alpha interi.
Berger-Parker è un altro booleano che permette di calcolare il corrispondente indice.
È inoltre possibile cambiare la base del logaritmo per il calcolo dell'indice di Rényi.

Spero di averti chiarito la struttura del codice.
Ti ringrazio per la grande disponibilità.
Martina
