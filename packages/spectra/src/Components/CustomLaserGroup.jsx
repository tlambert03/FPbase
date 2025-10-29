import { useApolloClient, useMutation, useQuery } from "@apollo/client"
import AddIcon from "@mui/icons-material/Add"
import DeleteIcon from "@mui/icons-material/Delete"
import Box from "@mui/material/Box"
import Button from "@mui/material/Button"
import IconButton from "@mui/material/IconButton"
import React, { useEffect, useRef, useState } from "react"
import { GET_EX_NORM, SET_EX_NORM, UPDATE_ACTIVE_SPECTRA } from "../client/queries"
import CustomLaserCreator from "./CustomLaserCreator"
import { categoryIcon } from "./FaIcon"

const CustomLaserGroup = React.memo(function CustomLaserGroup({ activeSpectra }) {
  const laserCounter = useRef(0)
  const [customLasers, setLasers] = useState([])
  const [updateSpectra] = useMutation(UPDATE_ACTIVE_SPECTRA)
  const {
    data: {
      exNorm: [, normID],
    },
  } = useQuery(GET_EX_NORM)

  const client = useApolloClient()
  const setExNorm = React.useCallback(
    (data) => client.mutate({ mutation: SET_EX_NORM, variables: { data } }),
    [client]
  )

  const clearNorm = React.useCallback(
    () =>
      client.mutate({
        mutation: SET_EX_NORM,
        variables: { data: [null, null] },
      }),
    [client]
  )

  useEffect(() => {
    if (activeSpectra && activeSpectra.length > 0) {
      const newLasers = activeSpectra.filter(
        (as) =>
          as.startsWith("$cl") && !customLasers.find((item) => item.startsWith(as.split("_")[0]))
      )
      if (newLasers.length) {
        const inds = newLasers.map((id) => +id.split("_")[0].replace("$cl", ""))
        laserCounter.current = Math.max(...inds) + 1
        setLasers([...customLasers, ...newLasers])
      }
    }
  }, [activeSpectra, customLasers]) // eslint-disable-line

  const addRow = () => {
    setLasers([...customLasers, `$cl${laserCounter.current++}`])
  }

  const removeRow = (laser) => {
    const laserID = laser.split("_")[0]
    if (laserID === normID) {
      clearNorm()
    }
    setLasers(customLasers.filter((id) => !id.startsWith(laserID)))
    updateSpectra({
      variables: {
        remove: [laserID],
      },
    })
  }

  return (
    <div>
      {customLasers.sort().map((laser) => (
        <div style={{ width: "100%", margin: "4px 0" }} key={laser}>
          <Box display="flex" alignItems="center">
            {categoryIcon("CL", "rgba(0,0,50,0.4)", {
              style: {
                position: "relative",
                top: 0,
                left: 2,
                height: "1.3rem",
                marginRight: 10,
              },
            })}
            <Box flexGrow={1}>
              <CustomLaserCreator
                key={laser.split("_")[0]}
                id={laser}
                setExNorm={setExNorm}
                clearNorm={clearNorm}
                normID={normID}
              />
            </Box>
            <Box>
              <IconButton
                aria-label="Delete"
                color="secondary"
                tabIndex={-1}
                onClick={() => removeRow(laser)}
                style={{
                  padding: "6px 6px",
                  marginRight: 2,
                  marginLeft: 2,
                }}
              >
                <DeleteIcon />
              </IconButton>
            </Box>
          </Box>
        </div>
      ))}
      <Button
        variant="contained"
        color="primary"
        onClick={() => addRow()}
        style={{ marginTop: 8, marginLeft: 34 }}
      >
        <AddIcon />
        {`Add Laser`}
      </Button>
    </div>
  )
})

export default CustomLaserGroup
