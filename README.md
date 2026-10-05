# kdb-proj-hft
High-Frequency Market Microstructure &amp; Orderbook Engine (L2/L3) in KDB+/q
Un projet end-to-end de simulation de marché et d'analytics de microstructure temporelle en temps réel et historique

Architecture du projet (Sur ton GitHub)Feedhandler C++/Python (PyKX) :

Un simulateur/rejoueur de flux L2/L3 (p. ex. rejouer des données ITCH/OUCH ou L2 de Binance/Crypto via WebSockets ou de données historiques NASDAQ).   

Real-Time Orderbook Engine (.q) :Reconstruction in-memory du carnet d'ordres jusqu'au niveau 10 (Depth of Book) à partir des événements de type Add, Cancel, Execute.Maintien dynamique de l'état du carnet et calcul d'indicateurs de microstructure à haute fréquence à chaque tick :Orderbook Imbalance (déséquilibre bid/ask pondéré par la profondeur).Micro-price & Mid-price Drift.Effective & Realized Spread.

Stream Analytics & Feature Store (RTE - Real-Time Engine) :Agrégation temps réel avec aj / asof joins et fenêtres glissantes (xbar).   Moteur de détection d'anomalies / signaux : détection de Spoofing / Large Trades ou de Toxic Flow (VPIN - Volume-Synchronized Probability of Toxicity).

Storage & Optimisation HDB :Partitionnement par date/heure, compression ZSTD, utilisation stratégique des attributs (`p#, `s#).   Script de post-traitement EOD (End Of Day) sous PyKX/q.   

API & Dashboard (Python / PyKX Gateway) :Une Gateway Python/FastAPI interrogeant kdb+ via PyKX pour exposer des endpoints REST / WebSockets vers un dashboard léger (Streamlit ou Plotly) affichant le carnet d'ordres animé et les métriques de risque/PnL. 
