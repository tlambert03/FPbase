[![Logo](fpbase/static/images/logo_green_wide@1x.png)](https://www.fpbase.org)

# FPbase: The Fluorescent Protein Database

[![License: MIT](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![cookie](https://img.shields.io/badge/built%20with-Cookiecutter%20Django-brightgreen.svg)](https://github.com/pydanny/cookiecutter-django/)

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

#### TODO

* mutations / lineages
* FRET views and FRET partners
* weighted scoring based on user-selected attributes
* user scoring/comments
* user submit unpublished data
* 2P absorption cross-section
* photoswitchable chart
* bleach comparisons within studies
* scale spectra to EC/QY/brightness/excitation line
* allow spectra/attributes to be embedded elsewhere
* "add to comparison" option, with comparison page, sequence alignment
* molecular weight
* OSER data
* Maturation data
