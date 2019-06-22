import React from "react"
import IconButton from "@material-ui/core/IconButton"
import SettingsIcon from "@material-ui/icons/Settings"
import SwipeableDrawer from "@material-ui/core/SwipeableDrawer"
import OwnerOptionsForm from "./OwnerOptionsForm"
import ChartOptionsForm from "./SpectraViewer/ChartOptionsForm"
import Button from "@material-ui/core/Button"
import CloseIcon from "@material-ui/icons/Close"
import DeleteIcon from "@material-ui/icons/Cached"
import { makeStyles } from "@material-ui/core/styles"

export const useStyles = makeStyles(theme => ({
  root: {
    width: "auto",
    padding: 30,
    paddingBottom: 20,
    [theme.breakpoints.down("sm")]: {
      paddingTop: 20
    }
  }
}))

const SettingsDrawer = ({ clearForm }) => {
  const classes = useStyles()
  const [drawerOpen, setDrawerOpen] = React.useState(false)

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
            <CloseIcon></CloseIcon>
          </Button>
          <Button
            variant="contained"
            color="secondary"
            style={{ float: "right" }}
            onClick={() => {
              clearForm()
              setDrawerOpen(false)
            }}
          >
            <DeleteIcon /> &nbsp; Remove All Spectra
          </Button>
        </div>
      </SwipeableDrawer>
    </div>
  )
}

export default SettingsDrawer
