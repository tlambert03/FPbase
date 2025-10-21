import "./js/detect-touch"
import { createElement } from "react"
import { createRoot } from "react-dom/client"
import App from "@fpbase/spectra"

const root = createRoot(document.getElementById("spectra-viewer"))
root.render(createElement(App, { uri: "/graphql/" }, null))
