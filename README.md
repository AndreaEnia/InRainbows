# In Rainbows: Spatially-resolved SED fitting of local galaxies in Dustpedia

## Introduction

To run a galaxy (i.e. NGCxxxx), there are three simple steps to follow:

- Create the galaxy folder, with the same name as the galaxy (NGCxxxx, mandatory)
- cut'n'paste inside the folder the notebook "GalaxyFlow.ipynb" found in scripts/
- follow the instructions inside the notebook

The notebook is a series of simple istructions, as it creates step by step the necessary folders with the reduced maps, results and plots. Check under for an exhaustive expalanation of what each function does, and how to customize them, in the meantime a general overview.

**Preliminary steps.**

1) Download from the DustPedia website bands and associated errors (if available) the user wants to work with. I suggest to pre-check in the site the bands and errors availability.

2) Data reduction, and for data reduction I mean: sky subtraction (if present), degradation to the worst map resolution, which as of today is assumed to always be SPIRE350. Yes, this thing could easily be generalized, but I don't tengo cazzi. This step is time consuming, so the reduced maps are saved in a folder, as in this way there is no need to always redo the whole procedure.

3)  Measure the maps rms (per pixel).

**Photometric steps.**

4) The user, depending on what results it desires, has to simply fill there four variables:
```python
* run_type = # The run name, i.e. pBp for pixel-by-pixel, or 1kpc for 1kpc apertures, o "Bh<oalf3-2Hhlad" if it feels creative.
* aperture_side = # How big the SIDE of the (square) apertures should be, either in parsec, kiloparsec or arcsec.
* reference_map_4grid = # The reference map for the apertures grid generation (coordinates.txt with values in ra e dec)
* reference_map_4check = # The reference map to check how the apertures were created 
* grid_extention_from_center = # How long the apertures grid should be extended from the galaxy center, either in pc, kpc or arcsec
* plot_radius = # How big the grid check plot should be. I suggest slightly bigger than the former grid_extention
```

5) Now it is the time for the big cell with all the photometric steps: the grid is generated, checked (eventually) A questo punto c'è il megacellone con gli step fotometrici: si genera la griglia, la si controlla (eventualmente), si fa la fotometria su quella griglia di aperture, si costruisce la tabella che andrà dentro magphys per il SED fitting. Quest'ultima non viene costruita su tutte le aperture, MA su quelle che superano una treshold di punti fotometrici positivi (di default 10), considerato che se il SNR va sotto un'altra threshold (di default 1.5), quel punto è escluso dalla fotometria.

**Step di SED fitting e risultati.**

6) A questo punto c'è solo da copiare la tabella dentro la cartella di magphys, farlo girare, aspettare che finisca, ricopiare i risultati dentro la cartella 'magphys'+_run_type e leggerli, salvando i risultati in un comodo magphys_results.csv

7) Uno step aggiuntivo opzionale, vista l'eventuale presenza di punti sparsi sconnessi dal resto delle aperture per via della treshold che riportavo prima. Questo step è il clustering removal, una roba ipermega tecnologica che cerca di salvare solo i punti che appartengono al grappolo principale di aperture, scartando quelli troppo isolati dal resto.

8) Infine, plotta i risultati. Questi plot sono molto brutti e preliminari, giusto per vedere cosa usciva fuori.

Descrizione pedante dei singoli scripts.
-----------
Tutto il lavoro descritto sopra viene eseguito in diversi scripts, tutti ospitati dentro la cartella scripts/ (duh!), i cui titoli sono abbastanza self-explanatory:
* `GenericUsefulScripts` (GUS)
* `DataReduction` (DataRed)
* `PhotometryScripts` (PhotoScripts)
* `SedFittingScripts` (SEDfit)
* `ResultsScripts` (ResScripts)

Iniziamo.
### GUS
GUS è un contenitore di un po' di tutto: gli script che scaricano le mappe (`GUS.download_data`, `GUS.download_data_errors`), le classi per leggere le proprietà fisiche delle galassie (`GUS.GalaxyProperties`), per leggere le mappe (`GUS.FitsUtils`) e per convertire le scale fisiche in quelle angolari e viceversa (`GUS.AngularPhysicalConv`), script per evitare ambiguità di nome per le bande osservative (`GUS.band_disambiguation`) e appioppare una colomap alle bande (`GUS.associate_colormap`), robe così va.

### DataRed
DataRed è lo script dei processi di riduzione dati e band degradation (e per verificare cosa è uscito fuori).
Per ridurli, si chiama così: `DataRed.data_reduction(galaxy_name, path_fits_input = 'Caapr/Maps')`. Il path di input dipende dal se si voglia o meno far girare CAAPR per eseguire la sottrazione delle stelle, come riportato sopra. Siccome non lo faccio girare mai, in automatico ho messo il path in cui vengono scaricate le mappe. Questo script genera una cartella *_ReducedMaps* nel quale vengono salvate le mappe ridotte.

Il check viene fatto con `DataRed.check_Dustpedia(galaxy_name, working_bands)`, e banalmente prende il flusso delle mappe ridotte entro l'apertura usata dalla collaborazione di Dustpedia, e lo confronta con il flusso riportato da loro (tutte info che stanno dentro la cartella *DustPedia_Tables*)

Nel caso in cui una mappa ridotta non vada bene (flussi totalmente scazzati), o si senta il bisogno di far rigirare il processo di riduzione, l'utente non deve spaventarsi, lo script salta in automatico le bande per cui trova la mappa ridotta dentro la cartella. Nel caso, quindi, bisognerà solo cancellare la mappa che non va bene, e far rigirare la cella.

### PhotoScripts
PhotoScripts sono tutti quegli script per stimare la rms sulle varie mappe, generare la griglia di coordinate e farci fotometria sopra, quindi costruire le tabelle fotometriche di interesse, tra cui quella da dare in pasto a MAGPHYS.

Per quanto riguarda la **rms**, dentro lo stesso script `PhotoScripts.evaluate_rms(galaxy_properties)` ci sono due versioni, v1 e v2. Entrambe le versioni condividono il grosso dello script, solo la parte finale cambia. Per grosso intendo: il mascheramento della galassia, in base a quando presente nelle tabelle di DustPedia, il sigma clipping della rimanente mappa, e la generazione di 2500 aperture random di 10 arcsec di raggio (di default, ma si possono fornire valori, basta mettere tra gli aromenti N_ap = xxx e ap_radius_physical = xxx). A questo punto i metodi divergono: v2 è quella che ho usato per i nostri paper, v1 è quella che uso adesso che sto a Bologna:
- v1 invece prende la distribuzione del flusso di tutte le 2500 aperture, e misura la rms come deviazione standard di questa distribuzione
- v2 è sostanzialmente copincollata dal modo in cui Chris Clark stima la rms dentro CAAPR. Clark misura la rms in tutte le 2500 aperture, e poi prende come rms della mappa il picco della distribuzione
In entrambi i casi la rms è PER PIXEL, quindi per misurare l'errore su una misura di flusso bisognerà moltiplicarla per il numero di pixel dentro l'apertura.

A questo punto inizia la fase fotometrica vera e propria, in cui è possibile customizzare al massimo il processo. Si riempiono le caselle di cui sopra, e si fanno girare in serie:
- lo script che genera la griglia di coordinate `PhotoScripts.generate_coordinate_grid_within_DustAp`. Questa griglia viene generata (in un modo ancora molto artigianale che devo migliorare) dentro cerchio di raggio fornito da te dal centro della galassia, ed è possibile anche fornire un ellitticità e un angolo di posizione, nel caso in cui i.e. la galassia sia molto schiacciata. Within Dustpedia Aperture nel senso che se uno tra angolo, ellitticità e semiasse di avoidance sono messi come False, lui piglia i valori per la singola galassia forniti nelle tabelle di DustPedia, ovvero quelli entro cui loro hanno fatto la fotometria;
- siccome capita di farla troppo grossa, o troppo piccola, c'è lo script per verificare: `PhotoScripts.check_coordinate_grid`;
- infine, si fa la fotometria mappa per mappa sulle aperture date con `PhotoScripts.do_photometry_build_table`. E' possibile fornire due threshold: una di SNR (1.5 di default), per cui tutti i valori con rapporti segnale rumore inferiore alle threshold vengono scartati e piazzati come -99, una di numero di fotometrie negative oltre il quale l'intera apertura viene scartata dal processo di SED fitting (10 di default), cosa che succede spesso nelle parti esterne delle galassie. Questo script genera tre tabelle **photometries_table.csv**, con tutte le fotometrie misurate anche nei punti che non soddisfano la threshold, **GOOD_photometries_table.csv**, scartando le aperture che non soddisfano il criterio, **Final_for_MAGPHYS.csv**, la tabella formattata nel formato che piace a MAGPHYS.

### SEDfit
Adesso le cose si fanno **DAVVERO** artigianali. E tutto perché la versione primigenia di MAGPHYS, per quanto bella e funzionale, non è molto user friendly. O meglio, le è nel senso che vuole solo due tabelle in input, quella dei flussi e quella dei filtri, ma non lo è nel senso che gira su tcsh e non su bash, non ha un sistema per parallelizzare le computazioni, e gli output vengono sparati direttamente nella cartella di MAGPHYS senza poter fornire un path, ergo **NON si possono far girare più galassie contemporameante**. Di conseguenza, queste ultime cose dovrò sistemarle io, a mano.
Questi script qui, `SEDfit.prepare_magphys`, `SEDfit.run_magphys`, `SEDfit.move_results` fanno sostanzialmente quanto appena detto.

**!!! DEVO PERO PRECISARE DEI CAVEAT, PERCHE' A SBAGLIARE QUALCOSA QUI SI PUO' IMPALLARE IL CODICE, O L'INTERO ATACAMA !!!**

- `SEDfit.prepare_magphys` banalmente crea dentro /home/dustpedia/magphys una cartella con il nome della galassia, e una sottocatella con il nome della run, dentro cui piazza il file con i filtri di interesse **filters.dat**, e spezza la singola tabella **Final_for_MAGPHYS.csv** in tante sottotabelle **observations_xxx_xxx.dat** a seconda di quante run parallele vuole fare l'utente: fill_cores, almost_fill_cores (tutti i cores meno due), single (niente parallelizazione), oppure un numero deciso dall'utente di entries per file. Se quest'ultimo fa si che in teoria si occupino più cores di quanti disponibili, il codice ti dice puppa. Infine, se necessario, si possono escludere dei filtri dal SED fitting.
- `SEDfit.run_magphys` fa girare MAGPHYS. Ora, sembra una cagata, ma è la cosa che, non essendo avvezzo a programmare in bash, mi ha preso più tempo. Perché? Perché MAGPHYS non è pensato per girare su una singola galassia più volta, ma per girare su un campione di diverse galassie a diverso z. Quindi funziona con questa sequenza: genera la libreria di redshift del campione, genera le librerie star+dust per ogni redshift associato, fitta il sample con quelle librerie. Ma io qui ho anche 10k punti allo stesso z, e lui non capisce che è la stessa galassia, ergo tenderebbe a generarmi 10k librerie stelle+dust sovrascrivendomele ogni volta. La cosa è palesemente demenziale. Ergo, devo fare questo giro della madonna: mi creo un obs_zero.dat con solo una riga, su questo genero le librerie star+dust, **ASPETTO CHE LE GENERI** con un brutalissimo "fermati per 300 secondi, poi fai il resto", quindi apro uno screen per ogni run parallela, modifico il file di configurazione .magphys_tcshrc, e mando il comando fit_sample. Tutto questo viene fatto in automatico dallo script.
- `SEDfit.run_magphys` controlla che MAGPHYS abbia finito, e poi sposta i risultati sia nella cartella della galassia creata su /home/dustpedia/magphys, che nella cartella 'magphys_'+run_type dentro la cartella /home/dustpedia/DustFolder/NGCxxx/ (ti chiederà di confermare che ti piaccia il path), da cui poi verrano letti. Questa duplicazione mi serviva per pararmi il culo da cancellazioni improvvide e accidentali. Come fa a controllare che MAGPHYS abbia finito? Mentre gira, MAGPHYS genera dei file di dimensione 0 bytes, xxxx.fit, dentro cui andrà il modello di bestfit. Finché ci sono file di 0 bytes, MAGPYHS non ha finito, e il codice ti dirà di tornare più tardi. Appena non ci sono più, MAGPHYS ha finito, e sposterà i risultati.

### ResScripts
Gli scripts contemporaneamente più semplici e incasinati. Semplici, perché con `ResScripts.read_magphys(galaxy_properties, run_type, aperture_side)` il codice legge tutti i vari xxx.sed e xxx.fit dentro la cartella con gli output di MAGPHYS, e da quelli genera un **NGC0628_results_pBp.dat** contenente tutti i risultati finali. Incasinati, perché oltre agli output standard di MAGPHYS, ci sono quelli di nostri interesse, come la distanza dalla MS o la distanza dal centro della galassia. Per tacere della fase di plot dei risultati, il che porta a personalizzare molto il codice. Di conseguenza, io mi prendevo la tabella di output, e poi a parte generavo i plot.

Ci tengo però a precisare il funzionamento di uno script che è lasciato lì di default, che è `clustering_removal(galaxy_properties, run_type, eps_value = 0.1, s_size = 10)`. Siccome capita spesso che ci siano dei punti sparsi sfuggiti ai criteri di thresholding, assolutamente sconnessi dal resto della galassia, e che poi nel plot della MS si manifestano come outliers imbarazzanti, ci tengo a farli fuori in virtù della loro sconnessione dal tessuto del grosso dei punti. Quindi con sklearn cerco di capire come sono clusterati i miei punti, o meglio, genero diversi cluster in base alla metrica punto/punto, esplicitata dentro quel parametro eps, che di default è preso pari a 0.09. Solitamente oscilla tra 0.06 (cerca di definire più clusters dentro il mio set di punti) e 0.5 (cerca di definire meno clusters dentro il mio set di punti). Partendo dal valore di default, lui genera un plot, s_size being la dimensione dell'apertura, purtroppo va settata a mano che non esiste un modo per settarla in automatico. Nel plot ci sono i cluster identificati, a sinistra, e come uscirebbe fuori il plot della SFR a destra scartando tutti i punti che ***non*** appartengono al main cluster. A quel punto ti chiederà se ti piace quello che stai vedendo, o no. Se il valore di default non va bene, e lo vedi perché vengono scartati troppi punti o addirittura viene considerato come main cluster soltanto il centro della galassia, dai in input *n*, e poi fornisci un nuovo valore di eps, a seconda che tu voglia definire più o meno clusters. Il processo è iterativo, finché non ti va bene, a quel punto lui piazza un chi quadro di 100 a tutti i punti non appartenenti al main cluster e salva tutto dentro una tabella **NGC0628_results_pBp_ClRm.dat**. Se vuoi uscire fuori dal processo iterativo, basta scrivere *break* in input.
