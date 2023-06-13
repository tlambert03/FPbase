import { createElement } from "react"
import { render } from "react-dom"
import App from "@fpbase/blast"

render(createElement(App, null, null), document.getElementById("blast-app"))
