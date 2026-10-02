import { createContext, useContext } from "react"

/** @type {import("react").Context<{ object: import("highcharts").Axis, id: string } | null>} */
export const AxisContext = createContext(null)

/** The chart's x axis, for components rendered alongside the chart */
export const useAxis = () => useContext(AxisContext)
