Caagle - A functional genomic database on *Candida glabrata* and close species
==============================================================================

#Prerequisites

To run this project, you need:

* Python >= 2.7 (but not 3.x)

* several Python libraries:

    * [Flask microframework](http://flask.pocoo.org)
    * [PyRNA](https://github.com/fjossinet/RNA-Science-Toolbox)
    * [PyMongo](https://api.mongodb.org/python/current/)

* [MongoDB](https://www.mongodb.org/)

* [Seqmap](http://www-personal.umich.edu/~jianghui/seqmap/)

#Quickstart

* start the MongoDB server

* import some [Candida Genome Database](http://www.candidagenome.org/) content by running the script import_CGD.py in the scripts directory

* launch the website by running the script init.py

* open your browser and go to http://127.0.0.1:5000


