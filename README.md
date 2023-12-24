[![Logo](_resources/logo_green_wide@1x.gif)](https://www.fpbase.org)

# FPbase: The Fluorescent Protein Database

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CI](https://github.com/tlambert03/FPbase/actions/workflows/ci.yml/badge.svg)](https://github.com/tlambert03/FPbase/actions/workflows/ci.yml)
[![Cov](https://codecov.io/gh/tlambert03/FPbase/branch/main/graph/badge.svg)](https://codecov.io/gh/tlambert03/FPbase)
[![DOI](https://zenodo.org/badge/DOI/10.1038/s41592-019-0352-8.svg)](https://doi.org/10.1038/s41592-019-0352-8)

Source code for <https://www.fpbase.org>

Documentation and info on using the site: <https://help.fpbase.org>.

See also: [Using FPbase: The Fluorescent Protein
Database](https://pubmed.ncbi.nlm.nih.gov/36107335/) (2023) *Methods Mol Biol* .
2023;2564:1-45. doi: 10.1007/978-1-0716-2667-2_1

## Installation for local development

1. Clone repo and cd into directory

    ```bash
    git clone https://github.com/tlambert03/FPbase.git
    cd FPbase
    ```

1. Create/activate environment **using python 3.11** with pipenv/virtualenv/conda
1. Install python requirements for local development

    ```bash
    pip install -r backend/requirements/local.txt
    # note: on mac silicon, you might have difficulty compiling psycopg2
    # in which case you should pip install psycopg2-binary instead
    ```

1. Make sure that you have `postgres` installed.

   On macOS, with homebrew:

   ```sh
   brew install postgresql@15
   brew services start postgresql@15
   ```

1. Install [Node.js](https://nodejs.org/en/) & [pnpm](https://pnpm.js.org/en/) (homebrew: `brew install node`)

    ```bash
    npm i -g pnpm
    ```

1. Install frontend requirements

    ```bash
    pnpm install
    ```

1. Install a local postgreSQL database (for mac: [postgres.app](https://postgresapp.com/))
1. Create database, and apply migrations

    ```bash
    createdb fpbase
    python backend/manage.py migrate
    ```

1. start dev servers:

    ```bash
    npm run start
    python backend/manage.py runserver
    ```

### How to cite FPbase

If you have used FPbase in a publication, or are referencing an FPbase protein
collection or microscope in your methods, please cite the following paper:

Lambert, TJ (2019) FPbase: a community-editable fluorescent protein database.
*Nature Methods*. doi:
[10.1038/s41592-019-0352-8](https://doi.org/10.1038/s41592-019-0352-8)

### Contributing

If you would like to contribute to the website directly (for instance, to add a
feature or fix an error), please branch off of develop and submit a pull
request.

If you have data that you would like to contribute to the database, please do
*not* do that here.  All data can be submitted directly on the website:

[Submit a fluorescent protein](https://www.fpbase.org/submit/)

[Submit spectral information](https://www.fpbase.org/spectra/submit/)

### Thank you to these providers for supporting open source projects

<br/>

[<img src="backend/proteins/static/images/logo-algolia-nebula-blue-full.svg"
width="170">](https://www.algolia.com/)

[<img src="_resources/sentry-logo-black.svg"
width="200">](https://sentry.io/)

[<img src="_resources/gitbook_avatar-rectangle.png"
width="200">](https://www.gitbook.com/)
