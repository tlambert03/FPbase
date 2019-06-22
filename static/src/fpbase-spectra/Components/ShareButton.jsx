import React, { useEffect, useState } from "react"
import Menu from "@material-ui/core/Menu"
import MenuItem from "@material-ui/core/MenuItem"
import ListItemIcon from "@material-ui/core/ListItemIcon"
import ListItemText from "@material-ui/core/ListItemText"
import DownloadIcon from "@material-ui/icons/GetApp"
import ChartIcon from "@material-ui/icons/InsertChart"
import PrintIcon from "@material-ui/icons/Print"
import LinkIcon from "@material-ui/icons/Link"
import CloseIcon from "@material-ui/icons/Close"
import Divider from "@material-ui/core/Divider"
import ShareIcon from "@material-ui/icons/Share"
import { useQuery } from "react-apollo-hooks"
import { GET_ACTIVE_SPECTRA, GET_CHART_OPTIONS } from "../client/queries"
import stateToUrl from "./stateToUrl"
import TextField from "@material-ui/core/TextField"
import Tooltip from "@material-ui/core/Tooltip"
import InputAdornment from "@material-ui/core/InputAdornment"
import IconButton from "@material-ui/core/IconButton"
import Highcharts from "highcharts"
import Dialog from "@material-ui/core/Dialog"
import DialogTitle from "@material-ui/core/DialogTitle"
import DialogContent from "@material-ui/core/DialogContent"
import DialogActions from "@material-ui/core/DialogActions"
import { FAIcon } from "./FaIcon"
import { faEnvelope, faCopy } from "@fortawesome/free-solid-svg-icons"
import { faTwitter } from "@fortawesome/free-brands-svg-icons"
import ClipboardJS from "clipboard"
import { makeStyles } from "@material-ui/core/styles"
import Zoom from "@material-ui/core/Zoom"

window.twttr = (function(d, s, id) {
  var js,
    fjs = d.getElementsByTagName(s)[0],
    t = window.twttr || {}
  if (d.getElementById(id)) return t
  js = d.createElement(s)
  js.id = id
  js.src = "https://platform.twitter.com/widgets.js"
  fjs.parentNode.insertBefore(js, fjs)

  t._e = []
  t.ready = function(f) {
    t._e.push(f)
  }

  return t
})(document, "script", "twitter-wjs")

const useStyles = makeStyles(theme => ({
  textField: {
    flexBasis: 200,
    width: "98%",
    margin: theme.spacing(1)
  }
}))

function ShareLinkAlert({ open, setOpen }) {
  const [qString, setQString] = useState("")
  const classes = useStyles()

  const {
    loading: spectraLoading,
    data: { activeSpectra }
  } = useQuery(GET_ACTIVE_SPECTRA)
  const {
    loading: chartLoading,
    data: { chartOptions }
  } = useQuery(GET_CHART_OPTIONS)

  useEffect(() => {
    if (!spectraLoading && !chartLoading) {
      setQString(stateToUrl(activeSpectra, chartOptions))
    }
  }, [activeSpectra, spectraLoading, chartOptions, chartLoading])

  const [tooltipOpen, setTooltipOpen] = React.useState(false)

  function handleTooltipOpen() {
    setTooltipOpen(true)
    setTimeout(() => setTooltipOpen(false), 1200)
  }

  useEffect(() => {
    const cp = new ClipboardJS("#copy-button", {
      target: function(trigger) {
        return document.getElementById("qString-input")
      }
    })
    cp.on("success", handleTooltipOpen)
    return () => {
      cp.destroy()
    }
  }, [])

  return (
    <div>
      <Dialog
        open={open}
        onClose={() => setOpen(false)}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
        fullWidth
        maxWidth={"md"}
      >
        <DialogTitle id="alert-dialog-title">
          You can use this URL to recreate the current graph:
        </DialogTitle>
        <DialogContent>
          <TextField
            id="qString-input"
            className={classes.textField}
            variant="outlined"
            type={"text"}
            label="URL"
            value={qString}
            fullWidth
            readOnly
            InputProps={{
              style: { color: "#669" },
              endAdornment: (
                <InputAdornment position="end">
                  <Tooltip
                    title="Copied!"
                    open={tooltipOpen}
                    disableFocusListener
                    disableHoverListener
                    disableTouchListener
                    placement="top"
                    TransitionComponent={Zoom}
                    TransitionProps={{ timeout: { enter: 200, exit: 700 } }}
                  >
                    <IconButton
                      edge="end"
                      id="copy-button"
                      aria-label="Toggle password visibility"
                    >
                      <FAIcon icon={faCopy} style={{}}></FAIcon>
                    </IconButton>
                  </Tooltip>
                </InputAdornment>
              )
            }}
          />
        </DialogContent>
        <DialogActions>
          <IconButton
            color="primary"
            href={`mailto:?&subject=Spectra%20at%20FPbase&body=${qString.replace(
              /&/g,
              "%26"
            )}`}
          >
            <FAIcon icon={faEnvelope} style={{}}></FAIcon>
          </IconButton>
          <IconButton
            color="primary"
            href={`http://twitter.com/intent/tweet?url=${qString.replace(
              /&/g,
              "%26"
            )}&text=Check out these spectra at FPbase%0D%0A&via=FPbase`}
            title="Tweet"
            target="_blank"
          >
            <FAIcon icon={faTwitter} style={{}}></FAIcon>
          </IconButton>
          <IconButton onClick={() => setOpen(false)} color="primary" autoFocus>
            <CloseIcon></CloseIcon>
          </IconButton>
        </DialogActions>
      </Dialog>
    </div>
  )
}

const ShareButton = () => {
  const [anchorEl, setAnchorEl] = useState(null)
  const [shareLinkOpen, setShareLinkOpen] = useState(false)
  const chart = Highcharts.charts[0]

  function handleShareClick(event) {
    setAnchorEl(event.currentTarget)
  }

  function handleClose() {
    setAnchorEl(null)
  }

  function exportChart(format) {
    if (format === "csv") {
      chart.downloadCSV()
    } else {
      chart.exportChart({
        type: format,
        filename: "FPbaseSpectra"
      })
    }
    setAnchorEl(null)
  }

  function printChart(format) {
    chart.print()
    setAnchorEl(null)
  }

  function openShareModal() {
    setShareLinkOpen(true)
    setAnchorEl(null)
  }

  const hasSeries = chart && chart.series.length > 0
  return (
    <div>
      <IconButton
        edge="end"
        color="inherit"
        aria-controls="simple-menu"
        aria-haspopup="true"
        onClick={handleShareClick}
        tabIndex={-1}
      >
        <ShareIcon />
      </IconButton>
      <Menu
        id="simple-menu"
        anchorEl={anchorEl}
        keepMounted
        open={Boolean(anchorEl)}
        onClose={handleClose}
      >
        {hasSeries ? (
          <div>
            <MenuItem onClick={() => exportChart("image/svg+xml")}>
              <ListItemIcon>
                <DownloadIcon />
              </ListItemIcon>
              <ListItemText primary="Download chart as SVG vector graphic" />
            </MenuItem>
            <MenuItem onClick={() => exportChart("image/png")}>
              <ListItemIcon>
                <DownloadIcon />
              </ListItemIcon>
              <ListItemText primary="Download chart as PNG image" />
            </MenuItem>
            <MenuItem onClick={() => exportChart("application/pdf")}>
              <ListItemIcon>
                <DownloadIcon />
              </ListItemIcon>
              <ListItemText primary="Download chart as PDF document" />
            </MenuItem>
            <Divider />
            <MenuItem onClick={() => exportChart("csv")}>
              <ListItemIcon>
                <ChartIcon />
              </ListItemIcon>
              <ListItemText primary="Download data as CSV" />
            </MenuItem>
            <Divider light />
            <MenuItem onClick={printChart}>
              <ListItemIcon>
                <PrintIcon />
              </ListItemIcon>
              <ListItemText primary="Print Chart" />
            </MenuItem>
            <Divider />
            <MenuItem onClick={openShareModal}>
              <ListItemIcon>
                <LinkIcon />
              </ListItemIcon>
              <ListItemText primary="Share Chart as URL" />
            </MenuItem>
            <ShareLinkAlert
              open={shareLinkOpen}
              setOpen={setShareLinkOpen}
            ></ShareLinkAlert>
          </div>
        ) : (
          <MenuItem disabled>
            <ListItemText primary="Add data to chart to enable exporting" />
          </MenuItem>
        )}
      </Menu>
    </div>
  )
}

export default ShareButton
