import React, { useEffect, memo, useState } from "react"
import Box from "@material-ui/core/Box"

import { UPDATE_ACTIVE_SPECTRA } from "../client/queries"
import { useMutation, useQuery } from "react-apollo-hooks"
import Typography from "@material-ui/core/Typography"
import InputSlider from "./InputSlider"
import { makeStyles } from "@material-ui/core/styles"
import Checkbox from "@material-ui/core/Checkbox"
import FormControlLabel from "@material-ui/core/FormControlLabel"
import RadioButtonUnchecked from "@material-ui/icons/RadioButtonUnchecked"
import RadioButtonChecked from "@material-ui/icons/RadioButtonChecked"

const useStyles = makeStyles({
  label: {
    fontSize: "small",
    color: "#444"
  }
})

// $cl1
const CustomLaserCreator = memo(function CustomLaserCreator({
  id,
  normID,
  setExNorm,
  clearNorm
}) {
  const classes = useStyles()

  let [laserID, _wave] = id.split("_")
  const [wave, setWave] = useState(_wave || 488)

  const updateSpectra = useMutation(UPDATE_ACTIVE_SPECTRA)
  useEffect(() => {
    updateSpectra({
      variables: {
        add: [`${laserID}_${wave}`],
        remove: [laserID]
      }
    })
    if (laserID === normID) {
      setExNorm([String(wave), laserID])
    }
    // for some reason setExNorm and updateSpectra are changing causing infinite loop when
    // adding a second laser line
    // eslint-disable-next-line
  }, [laserID, normID, wave])

  const handleNormCheck = (e, checked) => {
    if (checked) {
      setExNorm([String(wave), laserID])
    } else {
      clearNorm()
    }
  }

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
        <Box>
          <FormControlLabel
            value="normalize"
            label={
              <Typography className={classes.label}>
                Norm em. to this
              </Typography>
            }
            labelPlacement="end"
            control={
              <Checkbox
                icon={<RadioButtonUnchecked fontSize="small" />}
                checkedIcon={<RadioButtonChecked fontSize="small" />}
                value="checkedI"
                color="default"
                checked={id.startsWith(normID)}
                onChange={handleNormCheck}
              />
            }
          />
        </Box>
      </Box>
    </div>
  )
})

export default CustomLaserCreator
