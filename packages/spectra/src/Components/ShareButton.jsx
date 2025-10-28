import { useQuery } from '@apollo/client'
import { faTwitter } from '@fortawesome/free-brands-svg-icons'
import { faCopy, faEnvelope } from '@fortawesome/free-solid-svg-icons'
import CloseIcon from '@mui/icons-material/Close'
import DownloadIcon from '@mui/icons-material/GetApp'
import ChartIcon from '@mui/icons-material/InsertChart'
import LinkIcon from '@mui/icons-material/Link'
import PrintIcon from '@mui/icons-material/Print'
import ShareIcon from '@mui/icons-material/Share'
import Dialog from '@mui/material/Dialog'
import DialogActions from '@mui/material/DialogActions'
import DialogContent from '@mui/material/DialogContent'
import DialogTitle from '@mui/material/DialogTitle'
import Divider from '@mui/material/Divider'
import IconButton from '@mui/material/IconButton'
import InputAdornment from '@mui/material/InputAdornment'
import ListItemIcon from '@mui/material/ListItemIcon'
import ListItemText from '@mui/material/ListItemText'
import Menu from '@mui/material/Menu'
import MenuItem from '@mui/material/MenuItem'
import TextField from '@mui/material/TextField'
import Tooltip from '@mui/material/Tooltip'
import Zoom from '@mui/material/Zoom'
import { makeStyles } from '@mui/styles'
import ClipboardJS from 'clipboard'
import Highcharts from 'highcharts'
import React, { useEffect, useState } from 'react'
import { GET_ACTIVE_SPECTRA, GET_CHART_OPTIONS, GET_EX_NORM } from '../client/queries'
import { FAIcon } from './FaIcon'
import stateToUrl from './stateToUrl'

const useStyles = makeStyles((theme) => ({
  textField: {
    flexBasis: 200,
    width: '98%',
    margin: theme.spacing(1),
  },
  listIcon: {
    minWidth: 42,
  },
}))

function ShareLinkAlert({ open, setOpen }) {
  const [qString, setQString] = useState('')
  const classes = useStyles()

  const { loading: spectraLoading, data: spectraData } = useQuery(GET_ACTIVE_SPECTRA)
  const { loading: chartLoading, data: chartData } = useQuery(GET_CHART_OPTIONS)
  const { loading: exNormLoading, data: exNormData } = useQuery(GET_EX_NORM)

  const activeSpectra = spectraData?.activeSpectra || []
  const chartOptions = chartData?.chartOptions
  const exNorm = exNormData?.exNorm

  useEffect(() => {
    if (!spectraLoading && !chartLoading && !exNormLoading) {
      setQString(stateToUrl(activeSpectra, chartOptions, exNorm))
    }
  }, [activeSpectra, spectraLoading, chartOptions, chartLoading, exNormLoading, exNorm])

  const [tooltipOpen, setTooltipOpen] = React.useState(false)

  function handleTooltipOpen() {
    setTooltipOpen(true)
    setTimeout(() => setTooltipOpen(false), 1200)
  }

  useEffect(() => {
    const cp = new ClipboardJS('#copy-button', {
      target: (_trigger) => document.getElementById('qString-input'),
    })
    cp.on('success', handleTooltipOpen)
    return () => {
      cp.destroy()
    }
  }, [handleTooltipOpen])

  return (
    <div>
      <Dialog
        open={open}
        onClose={() => setOpen(false)}
        aria-labelledby="alert-dialog-title"
        aria-describedby="alert-dialog-description"
        fullWidth
        maxWidth="md"
      >
        <DialogTitle id="alert-dialog-title">
          You can use this URL to recreate the current graph:
        </DialogTitle>
        <DialogContent>
          <TextField
            id="qString-input"
            className={classes.textField}
            variant="outlined"
            type="text"
            label="URL"
            value={qString}
            fullWidth
            readOnly
            InputProps={{
              style: { color: '#669' },
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
                    <IconButton edge="end" id="copy-button" aria-label="Toggle password visibility">
                      <FAIcon icon={faCopy} style={{}} />
                    </IconButton>
                  </Tooltip>
                </InputAdornment>
              ),
            }}
          />
        </DialogContent>
        <DialogActions>
          <IconButton
            color="primary"
            href={`mailto:?&subject=Spectra%20at%20FPbase&body=${qString.replace(/&/g, '%26')}`}
          >
            <FAIcon icon={faEnvelope} style={{}} />
          </IconButton>
          <IconButton
            color="primary"
            href={`http://twitter.com/intent/tweet?url=${qString.replace(
              /&/g,
              '%26'
            )}&text=Check out these spectra at FPbase%0D%0A&via=FPbase`}
            title="Tweet"
            target="_blank"
          >
            <FAIcon icon={faTwitter} style={{}} />
          </IconButton>
          <IconButton onClick={() => setOpen(false)} color="primary" autoFocus>
            <CloseIcon />
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
  const classes = useStyles()

  function handleShareClick(event) {
    setAnchorEl(event.currentTarget)
  }

  function handleClose() {
    setAnchorEl(null)
  }

  function handleExportError(err, context) {
    console.error('Chart export failed:', err)
    if (window.Sentry) {
      window.Sentry.captureException(err, {
        tags: { component: 'SpectraViewer' },
        extra: { context, format: context },
      })
    }
    alert('Export failed. Please try again or contact us if the problem persists.')
  }

  function exportChart(format) {
    if (format === 'csv') {
      chart.downloadCSV()
    } else if (format === 'image/svg+xml') {
      // Safari-specific workaround: get SVG synchronously and download
      // Safari blocks downloads that happen asynchronously, so we must
      // get the SVG and trigger download in the same call stack

      // The timing is critical: Safari allows operations for up to 1 second after
      //  a user gesture, but:
      // - chart.exportChart() returns a Promise
      // - The actual download happens after the Promise resolves
      // - This asynchronous delay breaks the connection to the original click event
      // - Safari silently blocks the download because it's no longer user-initiated
      try {
        const svg = chart.getSVG()
        const blob = new Blob([svg], { type: 'image/svg+xml;charset=utf-8' })
        const url = URL.createObjectURL(blob)

        const link = document.createElement('a')
        link.href = url
        link.download = 'FPbaseSpectra.svg'
        document.body.appendChild(link)
        link.click()
        document.body.removeChild(link)

        // Clean up the blob URL
        setTimeout(() => URL.revokeObjectURL(url), 100)
      } catch (err) {
        handleExportError(err, 'SVG export (Safari workaround)')
      }
    } else {
      chart.exportChart({
        type: format,
        filename: 'FPbaseSpectra',
      })
    }
    setAnchorEl(null)
  }

  function printChart() {
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
            <MenuItem onClick={() => exportChart('image/svg+xml')}>
              <ListItemIcon className={classes.listIcon}>
                <DownloadIcon />
              </ListItemIcon>
              <ListItemText primary="Download chart as SVG" style={{ paddingRight: 20 }} />
            </MenuItem>
            <MenuItem onClick={() => exportChart('image/png')}>
              <ListItemIcon className={classes.listIcon}>
                <DownloadIcon />
              </ListItemIcon>
              <ListItemText primary="Download chart as PNG" />
            </MenuItem>
            <MenuItem onClick={() => exportChart('application/pdf')}>
              <ListItemIcon className={classes.listIcon}>
                <DownloadIcon />
              </ListItemIcon>
              <ListItemText primary="Download chart as PDF" />
            </MenuItem>
            <Divider />
            <MenuItem onClick={() => exportChart('csv')}>
              <ListItemIcon className={classes.listIcon}>
                <ChartIcon />
              </ListItemIcon>
              <ListItemText primary="Download data as CSV" />
            </MenuItem>
            <Divider light />
            <MenuItem onClick={printChart}>
              <ListItemIcon className={classes.listIcon}>
                <PrintIcon />
              </ListItemIcon>
              <ListItemText primary="Print chart" />
            </MenuItem>
            <Divider />
            <MenuItem onClick={openShareModal}>
              <ListItemIcon className={classes.listIcon}>
                <LinkIcon />
              </ListItemIcon>
              <ListItemText primary="Share chart as URL" />
            </MenuItem>
            <ShareLinkAlert open={shareLinkOpen} setOpen={setShareLinkOpen} />
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
