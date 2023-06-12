import React, { useEffect, useRef } from "react"
import { useMutation, useQuery } from "@apollo/react-hooks"
import { provideAxis } from "react-jsx-highcharts"
import Input from "@material-ui/core/Input"
import gql from "graphql-tag"
import { Tooltip } from "@material-ui/core"
import { withStyles } from "@material-ui/core/styles"

const LightTooltip = withStyles(theme => ({
  tooltip: {
    backgroundColor: theme.palette.common.white,
    color: "rgba(0, 0, 0, 0.87)",
    boxShadow: theme.shadows[1],
    fontSize: 12,
    margin: "0 13px 7px",
  },
}))(Tooltip)

const CLASSES = {
  minInput: {
    fontWeight: "bold",
    fontSize: "0.75rem",
    width: 30,
    color: "#444",
    position: "absolute",
  },
  maxInput: {
    position: "absolute",
    fontWeight: "bold",
    fontSize: "0.75rem",
    width: 30,
    color: "#444",
  },
}

const MUTATE_CHART_EXTREMES = gql`
  mutation SetChartExtremes($extremes: [Float]!) {
    setChartExtremes(extremes: $extremes) @client
  }
`

const GET_CHART_EXTREMES = gql`
  {
    chartOptions @client {
      extremes
    }
  }
`

let counter = 0
const XRangePickers = ({ getAxis, getHighcharts, visible }) => {
  const {
    data: {
      chartOptions: {
        extremes: [min, max],
      },
    },
  } = useQuery(GET_CHART_EXTREMES)
  const [mutateExtremes] = useMutation(MUTATE_CHART_EXTREMES)
  const minNode = useRef()
  const maxNode = useRef()
  const axis = getAxis()
  const forceUpdate = React.useState()[1]

  useEffect(() => {
    if (min || max) {
      axis.setExtremes(min && Math.round(min), max && Math.round(max))
      if (min && max) {
        axis.object.chart.showResetZoom()
      }
    }
  }, []) // eslint-disable-line

  useEffect(() => {
    function handleAfterSetExtremes() {
      const e = axis.object.getExtremes()
      if (e) {
        const extremes = [
          e.userMin && Math.round(e.min),
          e.userMax && Math.round(e.max),
        ]
        // this seems to be causing a bug with the inputs
        mutateExtremes({ variables: { extremes } })
        forceUpdate(counter++)
      }
    }

    function positionInputs() {
      if (axis.object.labelGroup && minNode.current) {
        let leftPad = -5
        if (axis.object.chart.get("yAx1")) {
          leftPad += +axis.object.chart.get("yAx1").axisTitleMargin
        }
        let rightPad = 0
        if (axis.object.chart.get("yAx2")) {
          rightPad += +axis.object.chart.get("yAx2").axisTitleMargin
        }
        minNode.current.parentElement.style.left = `${leftPad}px`
        maxNode.current.parentElement.style.right = `${rightPad}px`
        axis.object.labelGroup.element.childNodes.forEach(
          node => (node.style.display = "block")
        )
        const { min: exMin, max: exMax } = axis.getExtremes()
        axis.object.labelGroup.element.childNodes.forEach(node => {
          if (
            Math.min(
              Math.abs(node.textContent - exMin),
              Math.abs(node.textContent - exMax)
            ) <
            0.43 * axis.object.tickInterval
          ) {
            node.style.display = "none"
          }
        })
      }
    }

    const Highcharts = getHighcharts()
    Highcharts.addEvent(axis.object.chart, "redraw", positionInputs)
    Highcharts.addEvent(axis.object, "afterSetExtremes", handleAfterSetExtremes)

    Highcharts.addEvent(
      axis.object.chart,
      "redraw",
      handleAfterSetExtremes
      // debounce(handleAfterSetExtremes, 200)
    )
    handleAfterSetExtremes()
    positionInputs()
    return () => {
      Highcharts.removeEvent(
        getAxis().object,
        "afterSetExtremes",
        handleAfterSetExtremes
        // debounce(handleAfterSetExtremes, 200)
      )
    }
  }, []) // eslint-disable-line

  const updateRange = () => {
    const extremes = [
      +minNode.current.value || null,
      +maxNode.current.value || null,
    ]
    axis.setExtremes(...extremes)
  }
  const handleKeyPress = e => {
    if (e.key === "Enter") {
      updateRange()
      e.target.select()
    } else {
      if (e.target.name === "min") {
        minNode.current.value = e.target.value
      }
      if (e.target.name === "max") {
        maxNode.current.value = e.target.value
      }
    }
  }

  const extremes = axis.getExtremes()
  const minColor = !extremes.userMin
    ? "444"
    : extremes.dataMin < extremes.min
    ? "#B1191E"
    : "#5F67CE"
  const maxColor = !extremes.userMax
    ? "444"
    : extremes.dataMax > extremes.max
    ? "#B1191E"
    : "#5F67CE"

  return (
    <div
      className="x-range-pickers"
      style={{
        height: 0,
        position: "relative",
        bottom: 38,
        display: visible ? "block" : "none",
      }}
    >
      <LightTooltip
        title="Type to change min, clear to autoscale"
        placement="top-end"
        TransitionProps={{ timeout: { enter: 150, exit: 400 } }}
      >
        <Input
          name="min"
          type="text"
          placeholder={`${extremes.dataMin || ""}`}
          value={Math.round(min) || ""}
          inputRef={minNode}
          onChange={e =>
            mutateExtremes({ variables: { extremes: [e.target.value, max] } })
          }
          onKeyPress={handleKeyPress}
          onBlur={updateRange}
          style={{ ...CLASSES.minInput, color: minColor }}
          inputProps={{ style: { textAlign: "center" } }}
        />
      </LightTooltip>
      <LightTooltip
        title="Type to change max, clear to autoscale"
        placement="top-start"
        TransitionProps={{ timeout: { enter: 150, exit: 400 } }}
      >
        <Input
          name="max"
          type="text"
          placeholder={`${extremes.dataMax || ""}`}
          value={Math.round(max) || ""}
          inputRef={maxNode}
          onChange={e =>
            mutateExtremes({ variables: { extremes: [min, e.target.value] } })
          }
          onKeyPress={handleKeyPress}
          onBlur={updateRange}
          style={{ ...CLASSES.maxInput, color: maxColor }}
          inputProps={{ style: { textAlign: "center" } }}
        />
      </LightTooltip>
    </div>
  )
}
export default provideAxis(XRangePickers)
