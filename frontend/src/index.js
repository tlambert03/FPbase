import "regenerator-runtime/runtime"
import "select2/dist/css/select2.css"
import "select2-theme-bootstrap4/dist/select2-bootstrap.css"
import "nouislider/distribute/nouislider.min.css"
import "./css/style.scss"
import "./css/nv.d3.css"

import "bootstrap"

import "select2/dist/js/select2.full.js"
import "./js/my-fontawesome.js"

import "./js/project.js"
import "./js/search_logic.js"
import "./js/protein_page.js"
import "./js/favit.js"
import "./js/jquery.formset.js"
import "./js/onload.js"
import "./js/microscope.js"
import "./js/scope_report.js"

import FPPropChart from "./js/ichart.js"
import initAutocomplete from "./js/algolia.js"
import LineageChart from "./js/lineage.js"
import initFRET from "./js/fret.js"
import * as Sentry from "@sentry/browser";

window.FPBASE = window.FPBASE || {}


(async () => {
  if (process.env.NODE_ENV === "production" && Boolean(process.env.SENTRY_DSN)) {
    try {
      await import("@sentry/browser");

      Sentry.init({
        dsn: process.env.SENTRY_DSN,
        release: process.env.SOURCE_VERSION,
        environment: process.env.NODE_ENV,
        // Session Replay
        integrations: [Sentry.replayIntegration()],
        replaysSessionSampleRate: 0.1,
        replaysOnErrorSampleRate: 1.0,
      })

      if (window.FPBASE.user) {
        Sentry.setUser(window.FPBASE.user);
      }

      window.Sentry = Sentry;
    } catch (error) {
      console.error("Failed to initialize Sentry:", error);
    }
  }
})();

window.FPBASE = {
  ...window.FPBASE,
  initAutocomplete,
  FPPropChart,
  LineageChart,
  initFRET
}
