import React, { memo } from "react"
import List from "@material-ui/core/List"
import { useMutation, useQuery } from "react-apollo-hooks"
import gql from "graphql-tag"
import { GET_CHART_OPTIONS } from "../../client/queries"
import { ListCheckbox } from "../ListCheckbox"
import Grid from "@material-ui/core/Grid"

const toggleMut = param => gql`
mutation Toggle${param} {
  toggle${param} @client
}
`

const ChartOptionsForm = memo(function ChartOptionsForm({ options }) {
  const toggleY = useMutation(toggleMut("YAxis"))
  const toggleX = useMutation(toggleMut("XAxis"))
  const toggleGrid = useMutation(toggleMut("Grid"))
  const toggleScaleEC = useMutation(toggleMut("ScaleEC"))
  const toggleScaleQY = useMutation(toggleMut("ScaleQY"))
  const toggleShareTooltip = useMutation(toggleMut("ShareTooltip"))
  const toggleAreaFill = useMutation(toggleMut("AreaFill"))
  const {
    data: { chartOptions }
  } = useQuery(GET_CHART_OPTIONS)

  return (
    <List dense>
      <Grid container spacing={3}>
        <Grid item sm={12} md={6} style={{margin: 0, padding: 0}}>
          <ListCheckbox
            onCheckItem={toggleY}
            checked={chartOptions.showY}
            name="Show Y Axis"
          />
          <ListCheckbox
            onCheckItem={toggleX}
            checked={chartOptions.showX}
            name="Show X Axis"
          />
          <ListCheckbox
            onCheckItem={toggleGrid}
            checked={chartOptions.showGrid}
            name="Show Grid"
          />
          <ListCheckbox
            onCheckItem={toggleAreaFill}
            checked={chartOptions.areaFill}
            name="Fill area under curves"
          />
        </Grid>
        <Grid item sm={12} md={6} style={{margin: 0, padding: 0}}>
          <ListCheckbox
            onCheckItem={toggleScaleEC}
            checked={chartOptions.scaleEC}
            name="Scale Excitation to Extinction Coefficient"
          />
          <ListCheckbox
            onCheckItem={toggleScaleQY}
            checked={chartOptions.scaleQY}
            name="Scale Emission to Quantum Yield"
          />
          <ListCheckbox
            onCheckItem={toggleShareTooltip}
            checked={chartOptions.shareTooltip}
            name="Share Tooltip (on hover)"
          />
        </Grid>
      </Grid>
    </List>
  )
})
export default ChartOptionsForm
