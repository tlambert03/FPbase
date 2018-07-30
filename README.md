[![Logo](fpbase/static/images/logo_green_wide@1x.gif)](https://www.fpbase.org)

# FPbase: The Fluorescent Protein Database

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![cookie](https://img.shields.io/badge/built%20with-Cookiecutter%20Django-brightgreen.svg)](https://github.com/pydanny/cookiecutter-django/)
[![Build Status](https://travis-ci.org/tlambert03/FPbase.svg?branch=develop)](https://travis-ci.org/tlambert03/FPbase)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1244328.svg)](https://doi.org/10.5281/zenodo.1244328)


https://www.fpbase.org

## Installation for local development

1. Clone repo and cd into directory
```
    $ git clone https://github.com/tlambert03/FPbase.git
    $ cd FPbase
```
3. Create/activate environment with pipenv/virtualenv/conda (python 3 required)
4. Install python requirements for local development
```
    $ pip install -r requirements/local.txt
```
5. Install [Node.js](https://nodejs.org/en/) & npm  (homebrew: `brew install node`)
6. Install frontend requirements and [gulpjs](https://gulpjs.com/)
```
    $ npm install
    $ npm install gulp-cli -g
```
7. Install a local postgreSQL database (for mac: [postgres.app](https://postgresapp.com/))
8. Create database, and apply migrations
```
    $ createdb fpbase
    $ python manage.py migrate
```
9. If desired, load sample data (this may eventually break if the database schema changes enough).
```
    $ python manage.py loaddata fixtures/testdata.json.gz
```
10. Compile assets and start server:
```
    $ gulp
```


### Contributing

Please branch off of develop for any pull requests.

### Cite as 

Talley Lambert. tlambert03/FPbase (2018). doi:10.5281/zenodo.1244328

#### TODO

* mutations / lineages
* weighted scoring based on user-selected attributes
* user scoring/comments
* user submit unpublished data
* photoswitchable chart
* bleach comparisons within studies
* allow spectra/attributes to be embedded elsewhere
* "add to comparison" option, with comparison page, sequence alignment
* molecular weight
* OSER data
* Maturation data
