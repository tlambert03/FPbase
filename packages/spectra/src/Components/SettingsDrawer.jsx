import DeleteIcon from "@mui/icons-material/Cached"
import CloseIcon from "@mui/icons-material/Close"
import SettingsIcon from "@mui/icons-material/Settings"
import Button from "@mui/material/Button"
import IconButton from "@mui/material/IconButton"
import SwipeableDrawer from "@mui/material/SwipeableDrawer"
import { makeStyles } from "@mui/styles"
import React, { useEffect } from "react"
import { useSpectraStore } from "../store/spectraStore"
import OwnerOptionsForm from "./OwnerOptionsForm"
import ChartOptionsForm from "./SpectraViewer/ChartOptionsForm"

export const useStyles = makeStyles((theme) => ({
  root: {
    width: "auto",
    padding: 30,
    paddingBottom: 20,
    [theme.breakpoints.down("sm")]: {
      paddingTop: 20,
    },
  },
}))

const SettingsDrawer = () => {
  const classes = useStyles()
  const [drawerOpen, setDrawerOpen] = React.useState(false)
  const clearAllSpectra = useSpectraStore((state) => state.clearAllSpectra)

  useEffect(() => {
    const handleKeyDown = (event) => {
      // don't do anything if we're on an input
      if (document.activeElement && document.activeElement.tagName.toUpperCase() === "INPUT") {
        return
      }
      switch (event.code) {
        case "Comma":
          event.preventDefault()
          setDrawerOpen(!drawerOpen)
          break
        default:
          break
      }
    }

    document.addEventListener("keydown", handleKeyDown)
    return () => {
      document.removeEventListener("keydown", handleKeyDown)
    }
  }, [drawerOpen])

  const handleClearForm = () => {
    clearAllSpectra()
    setDrawerOpen(false)
  }
  return (
    <div>
      <IconButton
        tabIndex={-1}
        onClick={() => setDrawerOpen(true)}
        edge="start"
        color="inherit"
        aria-label="Open drawer"
      >
        <SettingsIcon />
      </IconButton>
      <SwipeableDrawer
        anchor="bottom"
        open={drawerOpen}
        onClose={() => setDrawerOpen(false)}
        onOpen={() => setDrawerOpen(true)}
        BackdropProps={{ style: { backgroundColor: "rgba(0,0,0,0.2)" } }}
      >
        <div className={classes.root} role="presentation">
          <ChartOptionsForm />
          <OwnerOptionsForm />
          <Button
            color="secondary"
            onClick={() => setDrawerOpen(false)}
            style={{ marginLeft: "-25px" }}
          >
            <CloseIcon />
          </Button>
          <Button
            variant="contained"
            color="error"
            style={{ float: "right" }}
            onClick={handleClearForm}
          >
            <DeleteIcon /> &nbsp; Remove All Spectra
          </Button>
        </div>
      </SwipeableDrawer>
    </div>
  )
}

export default SettingsDrawer
