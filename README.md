[![Logo](static/src/images/logo_green_wide@1x.gif)](https://www.fpbase.org)

# FPbase: The Fluorescent Protein Database

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![cookie](https://img.shields.io/badge/built%20with-Cookiecutter%20Django-brightgreen.svg)](https://github.com/pydanny/cookiecutter-django/)
[![Build Status](https://travis-ci.org/tlambert03/FPbase.svg?branch=develop)](https://travis-ci.org/tlambert03/FPbase)
[![DOI](https://zenodo.org/badge/DOI/10.1038/s41592-019-0352-8.svg)](https://doi.org/10.1038/s41592-019-0352-8)

https://www.fpbase.org

documentation and info on using the site: https://help.fpbase.org

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

5. Install [Node.js](https://nodejs.org/en/) & npm (homebrew: `brew install node`)
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

### How to cite FPbase

If you have used FPbase in a publication, or are referencing an FPbase protein collection or microscope in your methods, please cite the following paper:

Lambert, TJ (2019) FPbase: a community-editable fluorescent protein database. _Nature Methods_. doi: [10.1038/s41592-019-0352-8](https://doi.org/10.1038/s41592-019-0352-8)

### Contributing

If you would like to contribute to the website directly (for instance, to add a feature or fix an error), please branch off of develop and submit a pull request.

If you have data that you would like to contribute to the database, please do _not_ do that here. All data can be submitted directly on the website:

[Submit a fluorescent protein](https://www.fpbase.org/submit/)

[Submit spectral information](https://www.fpbase.org/spectra/submit/)

### Thank you to these providers for supporting open source projects!

<br/>

[<img src="static/src/images/logo-algolia-nebula-blue-full.svg" width="170">](https://www.algolia.com/)

[<img src="static/src/images/sentry-logo-black.svg" width="200">](https://sentry.io/)

[<img src="static/src/images/TravisCI-Full-Color.png" width="190">](https://travis-ci.org/)

[<img src="static/src/images/gitbook_avatar-rectangle.png" width="200">](https://www.gitbook.com/)

[<img src="static/src/images/Browserstack-logo@2x.png" width="250">](https://www.browserstack.com)
