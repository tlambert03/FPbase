import React, { useEffect } from "react"
import Box from "@material-ui/core/Box"

import { UPDATE_ACTIVE_SPECTRA } from "../client/queries"
import { useMutation } from "react-apollo-hooks"
import Typography from "@material-ui/core/Typography"
import InputSlider from "./InputSlider"
import { makeStyles } from "@material-ui/core/styles"

export const useStyles = makeStyles({
  label: {
    fontSize: "small",
    color: "#444"
  }
})

// $cl1
const CustomLaserCreator = React.memo(function CustomLaserCreator({ id }) {
  const classes = useStyles()

  let [laserID, _wave] = id.split("_")
  const [wave, setWave] = React.useState(_wave || 488)

  const updateSpectra = useMutation(UPDATE_ACTIVE_SPECTRA)

  useEffect(() => {
    updateSpectra({
      variables: {
        add: [`${laserID}_${wave}`],
        remove: [laserID]
      }
    })
  }, [updateSpectra, laserID, _wave, wave])

  return (
    <div
      style={{
        padding: "10px 10px 0px 10px",
        border: "1px solid #ccc",
        borderRadius: 4
      }}
    >
      <Box display="flex" flexWrap="wrap">
        <Typography style={{ margin: "8px 10px 3px" }}>Laser</Typography>
        <Box flexGrow={2}>
          <div style={{ margin: "0 12px", minWidth: 200 }}>
            <Typography className={classes.label}>Wavelength</Typography>
            <InputSlider
              value={wave}
              setValue={setWave}
              valueLabelDisplay="auto"
              aria-labelledby="range-slider"
            />
          </div>
        </Box>
      </Box>
    </div>
  )
})

export default CustomLaserCreator
