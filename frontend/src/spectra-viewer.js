import "./js/detect-touch"
import { createElement } from "react"
import { render } from "react-dom"
import App from "@fpbase/spectra"

render(
  createElement(App, { uri: "/graphql/" }, null),
  document.getElementById("spectra-viewer")
)
