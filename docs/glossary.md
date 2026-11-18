# Glossary of German terms and abbreviations

Many identifiers, file names and comments in WaSiM are German or abbreviations of German words. This list translates the most frequent ones. Entries marked (?) are inferred from usage, not confirmed by the authors.

## Source files

| File | German | Meaning |
| --- | --- | --- |
| `wasimuzr.cc` | WaSiM + **u**ngesättigte **Z**one **R**ichards | main program of the Richards version |
| `modapuzr.cc` | **Mod**ell**ap**plikation uzr | model application: initialization and time loop |
| `modefuzr.cc` | **Mod**ell**ef**initionen (?) | base class `modelobject`: grids, memory, file I/O, MPI helpers |
| `defltuzr.cc` | **Def**au**lt**s uzr (?) | default/shared initialization: zones, regions, tables, MPI decomposition |
| `teilmod.cc` | **Teilmod**ell | submodel base class, statistics |
| `meteoint.cc` | **Meteo**rologische **Int**erpolation | meteorological interpolation |
| `expokorr.cc` | **Expo**sitions**korr**ektur | radiation/temperature correction for slope and aspect |
| `etpetr.cc` | ETP / ETR | potential / real evapotranspiration |
| `interzep.cc` | **Interzep**tion | interception |
| `snowmod.cc` | Schneemodell | snow and glacier model |
| `unsatzon.cc` | **unsat**urated **zon**e | unsaturated zone (Richards equation) |
| `gwflow.cc` | groundwater flow | 2D groundwater model |
| `abfluss.cc` | Abfluss | runoff routing, lakes, reservoirs, abstractions |
| `surfrout.cc`, `siltingup.cc` | | surface routing, silting up |

## Routing description keywords

| Keyword | German | Meaning |
| --- | --- | --- |
| `TG` | **T**eil**g**ebiet | subbasin |
| `OL` | **O**ber**l**ieger | upstream tributary subbasin routed into a `TG` |
| `SUMTRIB` | | sum of tributaries, routed in one channel |
| `ZL` | **Z**u**l**eitung | inflow, external or internal |
| `AL` | **A**b**l**eitung | abstraction, external or internal |
| `SP` | **S**peicher | reservoir or lake storage |
| `AE`, `AErel` | **A**nteil**E**inzugsgebiet (?) | (relative) catchment area in km² |
| `kh`, `kv` | | storage coefficient, main channel (**h**aupt) / flood plain (**V**orland) [h] |
| `Bh`, `Bv`, `Th` | **B**reite, **T**iefe | width of main channel / flood plain, depth of main channel [m] |
| `Mh`, `Mv` | **M**anning-Strickler | roughness, main channel / flood plain |
| `I`, `L` | | slope, length of the channel |

## Grid file extensions in the setups

| Extension | German | Meaning |
| --- | --- | --- |
| `.dhm`, `.dhk` | **d**igitales **H**öhen**m**odell | digital elevation model |
| `.ezg` | **E**in**z**ugs**g**ebiete | subbasin (zone) grid |
| `.fzs` | **F**ließ**z**eit**s**umme (?) | flow travel time sums |
| `.lnk` | link | river links (channel network) |
| `.slp`, `.exp` | slope, **Exp**osition | slope angle, aspect |
| `.use`, `.soil` | | land use, soil types |
| `.kol` | **Kol**mation | clogging of the river bed (leakage) |
| `.aq1`, `.kx1`, `.ky1`, `.s01` | **Aq**uifer | aquifer thickness, conductivity x/y, storage coefficient, layer 1 |
| `.bh1`, `.bq1` | **B**oundary **h**/**q** | groundwater boundary: fixed head / fixed flux |
| `.gk1` | **G**rundwasser-**K**olmation (?) | leakage coefficient to layer 1 |
| `.wit`, `.dep` | width, depth | river width and depth |
| `.glc` | glacier | glacier cells |
| `.sd1`…`.sd4` | slide | gravitational snow slide fractions (N, W, E, S) |

## Frequent words in identifiers and comments

| German | English | Example identifiers |
| --- | --- | --- |
| Abfluss, Abflüsse | runoff, discharge | `abflussmodellklasse` |
| Anteil | fraction, share | `schneeanteil`, `Anteil_hell` (daylight fraction) |
| Anzahl, anz | number, count | `gridmaxanzahl`, `anz_ZL` |
| aus / ein | out / in | `nodataaus`, `ncolsein` |
| Boden, Bodenfeuchte | soil, soil moisture | |
| Datei, dateiname | file, file name | `steuerdatei` |
| Dicke | thickness | `dicken`, `a_dicke` |
| Einzugsgebiet | catchment | |
| Ende, Anfang | end, start | `programmende` |
| Entfernung, entf | distance | `maxentf` |
| Fehlwert | missing value (nodata) | `fehlwert` |
| Firn | firn | `FIRN1`, `qfirn` |
| Flussdichte | drainage density | `flussdichte` |
| Gesamt, Gesamtgebiet | total, whole domain | `Gesamtgrid` |
| Grenze, Grenz- | threshold, limit | `Tgrenz`, `theta_grenz` |
| Grundwasser, GW, GWSTAND | groundwater, groundwater level | `old_gwstand` |
| hell / dunkel | light (day) / dark (night) | `Anteil_dunkel`, `Rrest_day/night` |
| Höhe, Höhenmodell | elevation, DEM | `Hoehenmodell` |
| Inhalt | content | |
| InV, InVtag | **In**ter**v**all: time step length in minutes; InVtag: time steps per day | `InV/1440` |
| Jahr, Monat, Tag, Stunde | year, month, day, hour | `stundea`, `stundee` (suffix `a` = Anfang/start, `e` = Ende/end of the model run) |
| Klasse | class | `schneemodellklasse` |
| Mittel, Mittelwert | mean | `Teilgebietsmittelwert` (subbasin mean) |
| Niederschlag | precipitation | `_Niederschlag_` |
| Oberlieger | upstream neighbour | `Oberliegers` |
| Rückgefrieren | refreezing | |
| Schicht | layer | `startschicht`, `endschicht` |
| Schmelze | melt | |
| Schnee | snow | `schneemodell` |
| See | lake | |
| Sonne, Strahlung | sun, radiation | |
| Spalte, Zeile | column, row | |
| Speicher, speichererror | storage, memory error | |
| Spende | specific discharge (per area) | `SpendenBuffer` |
| Steuerdatei, steuername | control file, control file name | |
| Summe, Summengrid | sum, sum grid | `summenwritegrid` |
| Teilgebiet, Teilmodell | subbasin, submodel | `teilmodellklasse` |
| Verdunstung | evaporation, evapotranspiration | `verdunstungsmodell` |
| Wert, Werte | value(s) | `zonenwert`, `wertzufuegen` (add value) |
| Wichtung | weighting | |
| Wurzel, Wurzeltiefe, WURZ | root, root depth | |
| Zeit, Zeitintervall | time, time step | `maxzeit` |
| Zone, Zonen | zone (subbasin or statistics zone) | `zzone`, `zonmax`, `zonenid` |
| Zulauf | inflow | |
