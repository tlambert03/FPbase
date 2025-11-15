import Grid from "@mui/material/Grid"
import List from "@mui/material/List"
import MenuItem from "@mui/material/MenuItem"
import Select from "@mui/material/Select"
import { memo, useCallback } from "react"
import PALETTES from "../../palettes"
import { useSpectraStore } from "../../store/spectraStore"
import ListCheckbox from "../ListCheckbox"

const ChartOptionsForm = memo(function ChartOptionsForm() {
  const chartOptions = useSpectraStore((state) => state.chartOptions)
  const updateChartOptions = useSpectraStore((state) => state.updateChartOptions)

  const toggleY = useCallback(
    () => updateChartOptions({ showY: !chartOptions.showY }),
    [chartOptions.showY, updateChartOptions]
  )
  const toggleX = useCallback(
    () => updateChartOptions({ showX: !chartOptions.showX }),
    [chartOptions.showX, updateChartOptions]
  )
  const toggleGrid = useCallback(
    () => updateChartOptions({ showGrid: !chartOptions.showGrid }),
    [chartOptions.showGrid, updateChartOptions]
  )
  const toggleScaleEC = useCallback(
    () => updateChartOptions({ scaleEC: !chartOptions.scaleEC }),
    [chartOptions.scaleEC, updateChartOptions]
  )
  const toggleScaleQY = useCallback(
    () => updateChartOptions({ scaleQY: !chartOptions.scaleQY }),
    [chartOptions.scaleQY, updateChartOptions]
  )
  const toggleShareTooltip = useCallback(
    () => updateChartOptions({ shareTooltip: !chartOptions.shareTooltip }),
    [chartOptions.shareTooltip, updateChartOptions]
  )
  const toggleAreaFill = useCallback(
    () => updateChartOptions({ areaFill: !chartOptions.areaFill }),
    [chartOptions.areaFill, updateChartOptions]
  )
  const setPalette = useCallback((palette) => updateChartOptions({ palette }), [updateChartOptions])

  if (!chartOptions) return null

  return (
    <List dense className="settings-list">
      <Grid container spacing={3}>
        <Grid item sm={12} md={6} style={{ margin: 0, padding: 0 }}>
          <ListCheckbox
            onCheckItem={toggleY}
            checked={chartOptions.showY}
            name={
              <span>
                Show Y Axis
                <span className="mobile-hide small text-muted ms-3">(Y key)</span>
              </span>
            }
          />
          <ListCheckbox
            onCheckItem={toggleX}
            checked={chartOptions.showX}
            name={
              <span>
                Show X Axis
                <span className="mobile-hide small text-muted ms-3">(X key)</span>
              </span>
            }
          />
          <ListCheckbox
            onCheckItem={toggleGrid}
            checked={chartOptions.showGrid}
            name={
              <span>
                Show Grid
                <span className="mobile-hide small text-muted ms-3">(G key)</span>
              </span>
            }
          />
          <ListCheckbox
            onCheckItem={toggleAreaFill}
            checked={chartOptions.areaFill}
            name={
              <span>
                Fill area under curves
                <span className="mobile-hide small text-muted ms-3">(F key)</span>
              </span>
            }
          />
        </Grid>
        <Grid item sm={12} md={6} style={{ margin: 0, padding: 0 }}>
          <ListCheckbox
            onCheckItem={toggleScaleEC}
            checked={chartOptions.scaleEC}
            name={
              <span>
                Scale Excitation to Extinction Coefficient
                <span className="mobile-hide small text-muted ms-3">(E key)</span>
              </span>
            }
          />
          <ListCheckbox
            onCheckItem={toggleScaleQY}
            checked={chartOptions.scaleQY}
            name={
              <span>
                Scale Emission to Quantum Yield
                <span className="mobile-hide small text-muted ms-3">(Q key)</span>
              </span>
            }
          />
          <ListCheckbox
            onCheckItem={toggleShareTooltip}
            checked={chartOptions.shareTooltip}
            name={
              <span>
                Show Y value for all spectra on chart hover
                <span className="mobile-hide small text-muted ms-3">(T key)</span>
              </span>
            }
          />
          <div style={{ marginBottom: 20 }}>
            <span style={{ marginRight: 14, marginLeft: 14 }}>Color palette:</span>
            <Select
              value={chartOptions.palette}
              onChange={(e) => setPalette(e.target.value)}
              inputProps={{
                name: "color-palette",
                id: "color-palette-select",
              }}
            >
              {Object.keys(PALETTES).map((i) => (
                <MenuItem value={i} key={i}>
                  {PALETTES[i].name}
                </MenuItem>
              ))}
            </Select>
            <span className="mobile-hide small text-muted ms-3">(C key)</span>
          </div>
        </Grid>
      </Grid>
    </List>
  )
})
export default ChartOptionsForm
