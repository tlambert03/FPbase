import React, { memo } from "react"
import List from "@material-ui/core/List"
import { useMutation, useQuery } from "@apollo/react-hooks"
import gql from "graphql-tag"
import Grid from "@material-ui/core/Grid"
import MenuItem from "@material-ui/core/MenuItem"
import Select from "@material-ui/core/Select"
import { GET_CHART_OPTIONS } from "../../client/queries"
import ListCheckbox from "../ListCheckbox"
import PALETTES from "../../palettes"

const toggleMut = param => gql`
mutation Toggle${param} {
  toggle${param} @client
}
`
const MUTATE_PALETTE = gql`
  mutation SetPalette($palette: String!) {
    setPalette(palette: $palette) @client
  }
`

const ChartOptionsForm = memo(function ChartOptionsForm({ options }) {
  const [toggleY] = useMutation(toggleMut("YAxis"))
  const [toggleX] = useMutation(toggleMut("XAxis"))
  const [toggleGrid] = useMutation(toggleMut("Grid"))
  const [toggleScaleEC] = useMutation(toggleMut("ScaleEC"))
  const [toggleScaleQY] = useMutation(toggleMut("ScaleQY"))
  const [toggleShareTooltip] = useMutation(toggleMut("ShareTooltip"))
  const [toggleAreaFill] = useMutation(toggleMut("AreaFill"))
  const [mutatePalette] = useMutation(MUTATE_PALETTE)
  const {
    data: { chartOptions },
  } = useQuery(GET_CHART_OPTIONS)

  return (
    <List dense className="settings-list">
      <Grid container spacing={3}>
        <Grid item sm={12} md={6} style={{ margin: 0, padding: 0 }}>
          <ListCheckbox
            onCheckItem={toggleY}
            checked={chartOptions.showY}
            name={(
              <span>
                Show Y Axis
                <span className="mobile-hide small text-muted ml-3">
                  (Y key)
                </span>
              </span>
)}
          />
          <ListCheckbox
            onCheckItem={toggleX}
            checked={chartOptions.showX}
            name={(
              <span>
                Show X Axis
                <span className="mobile-hide small text-muted ml-3">
                  (X key)
                </span>
              </span>
)}
          />
          <ListCheckbox
            onCheckItem={toggleGrid}
            checked={chartOptions.showGrid}
            name={(
              <span>
                Show Grid
                <span className="mobile-hide small text-muted ml-3">
                  (G key)
                </span>
              </span>
)}
          />
          <ListCheckbox
            onCheckItem={toggleAreaFill}
            checked={chartOptions.areaFill}
            name={(
              <span>
                Fill area under curves
                <span className="mobile-hide small text-muted ml-3">
                  (F key)
                </span>
              </span>
)}
          />
        </Grid>
        <Grid item sm={12} md={6} style={{ margin: 0, padding: 0 }}>
          <ListCheckbox
            onCheckItem={toggleScaleEC}
            checked={chartOptions.scaleEC}
            name={(
              <span>
                Scale Excitation to Extinction Coefficient
                <span className="mobile-hide small text-muted ml-3">
                  (E key)
                </span>
              </span>
)}
          />
          <ListCheckbox
            onCheckItem={toggleScaleQY}
            checked={chartOptions.scaleQY}
            name={(
              <span>
                Scale Emission to Quantum Yield
                <span className="mobile-hide small text-muted ml-3">
                  (Q key)
                </span>
              </span>
)}
          />
          <ListCheckbox
            onCheckItem={toggleShareTooltip}
            checked={chartOptions.shareTooltip}
            name={(
              <span>
                Show Y value for all spectra on chart hover
                <span className="mobile-hide small text-muted ml-3">
                  (T key)
                </span>
              </span>
)}
          />
          <div style={{ marginBottom: 20 }}>
            <span style={{ marginRight: 14, marginLeft: 14 }}>
              Color palette:
            </span>
            <Select
              value={chartOptions.palette}
              onChange={e =>
                mutatePalette({ variables: { palette: e.target.value } })
              }
              inputProps={{
                name: "color-palette",
                id: "color-palette-select",
              }}
            >
              {Object.keys(PALETTES).map(i => (
                <MenuItem value={i} key={i}>
                  {PALETTES[i].name}
                </MenuItem>
              ))}
            </Select>
            <span className="mobile-hide small text-muted ml-3">(C key)</span>
          </div>
        </Grid>
      </Grid>
    </List>
  )
})
export default ChartOptionsForm
