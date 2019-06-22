import React, { useEffect, useRef } from "react"
import { useMutation, useQuery } from "react-apollo-hooks"
import { provideAxis } from "react-jsx-highcharts"
import Input from "@material-ui/core/Input"
import gql from "graphql-tag"

const CLASSES = {
  minInput: {
    fontWeight: "bold",
    fontSize: "0.75rem",
    width: 30,
    color: "#444",
    position: "absolute"
  },
  maxInput: {
    position: "absolute",
    fontWeight: "bold",
    fontSize: "0.75rem",
    width: 30,
    color: "#444"
  }
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
        extremes: [min, max]
      }
    }
  } = useQuery(GET_CHART_EXTREMES)
  const mutateExtremes = useMutation(MUTATE_CHART_EXTREMES)
  const minNode = useRef()
  const maxNode = useRef()
  const axis = getAxis()
  const forceUpdate = React.useState()[1];

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
          e.userMax && Math.round(e.max)
        ]
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
        const { min, max } = axis.getExtremes()
        axis.object.labelGroup.element.childNodes.forEach(node => {
          if (
            Math.min(
              Math.abs(node.textContent - min),
              Math.abs(node.textContent - max)
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
    Highcharts.addEvent(axis.object, 'afterSetExtremes', handleAfterSetExtremes);

    Highcharts.addEvent(
      axis.object.chart,
      "redraw",
      handleAfterSetExtremes
      //debounce(handleAfterSetExtremes, 200)
    )
    handleAfterSetExtremes()
    positionInputs()
    return () => {
      Highcharts.removeEvent(
        getAxis().object,
        "afterSetExtremes",
        handleAfterSetExtremes
        //debounce(handleAfterSetExtremes, 200)
      )
    }
  }, []) // eslint-disable-line

  const updateRange = () => {
    const extremes = [
      +minNode.current.value || null,
      +maxNode.current.value || null
    ]
    axis.setExtremes(...extremes)
  }
  const handleKeyPress = e => {
    if (e.key === "Enter") {
      updateRange()
      e.target.select()
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
        display: visible ? "block" : "none"
      }}
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
    </div>
  )
}
export default provideAxis(XRangePickers)
