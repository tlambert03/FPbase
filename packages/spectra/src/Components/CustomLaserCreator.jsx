import React, { useEffect, memo, useState } from "react"
import Box from "@material-ui/core/Box"

import { useApolloClient } from "@apollo/react-hooks"
import Typography from "@material-ui/core/Typography"
import { makeStyles } from "@material-ui/core/styles"
import Checkbox from "@material-ui/core/Checkbox"
import FormControlLabel from "@material-ui/core/FormControlLabel"
import RadioButtonUnchecked from "@material-ui/icons/RadioButtonUnchecked"
import RadioButtonChecked from "@material-ui/icons/RadioButtonChecked"
import InputSlider from "./InputSlider"
import { UPDATE_ACTIVE_SPECTRA } from "../client/queries"

const useStyles = makeStyles({
  label: {
    fontSize: "small",
    color: "#444",
  },
})

// $cl1
const CustomLaserCreator = memo(function CustomLaserCreator({
  id,
  normID,
  setExNorm,
  clearNorm,
}) {
  const classes = useStyles()

  const [laserID, _wave] = id.split("_")
  const [wave, setWave] = useState(_wave || 488)

  const client = useApolloClient()
  useEffect(() => {
    client.mutate({
      mutation: UPDATE_ACTIVE_SPECTRA,
      variables: {
        add: [`${laserID}_${wave}`],
        remove: [laserID],
      },
    })

    if (laserID === normID) {
      setExNorm([String(wave), laserID])
    }
  }, [client, laserID, normID, setExNorm, wave])

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
        borderRadius: 4,
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
            label={(
              <Typography className={classes.label}>
                Norm em. to this
              </Typography>
)}
            labelPlacement="end"
            control={(
              <Checkbox
                icon={<RadioButtonUnchecked fontSize="small" />}
                checkedIcon={<RadioButtonChecked fontSize="small" />}
                value="checkedI"
                color="default"
                checked={id.startsWith(normID)}
                onChange={handleNormCheck}
              />
)}
          />
        </Box>
      </Box>
    </div>
  )
})

export default CustomLaserCreator
