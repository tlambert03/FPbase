import React, { useCallback, useEffect, useMemo } from "react"
import MyAppBar from "./Components/MyAppBar"
import OwnersContainer from "./Components/OwnersContainer"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import useKeyboardShortcuts from "./Components/useKeyboardShortcuts"
import WelcomeModal from "./Components/WelcomeModal"
import { useSpectraMetadata } from "./hooks/useSpectraMetadata"
import { useOwnerInfo, useSpectraInfo } from "./store/metadataStore"
import { useSpectraStore } from "./store/spectraStore"
import { parseURLParams } from "./utils/urlParams"
import "./polyfills"

const daysSinceLaunch = Math.round(
  (new Date(2019, 6, 1).getTime() - Date.now()) / (1000 * 60 * 60 * 24)
)

const App = () => {
  // Parse URL params and apply to store on initial mount
  useEffect(() => {
    const urlState = parseURLParams(window.location.search)

    // Only apply if URL has params
    if (Object.keys(urlState).length > 0) {
      const store = useSpectraStore.getState()

      // Apply custom filters/lasers first (they define the stable IDs)
      if (urlState.customFilters) {
        for (const [id, params] of Object.entries(urlState.customFilters)) {
          store.addCustomFilter(id, params)
        }
      }
      if (urlState.customLasers) {
        for (const [id, params] of Object.entries(urlState.customLasers)) {
          store.addCustomLaser(id, params)
        }
      }

      // Then apply active spectra (now uses stable IDs like "$cf0", "$cl1")
      if (urlState.activeSpectra) {
        store.setActiveSpectra(urlState.activeSpectra)
      }
      if (urlState.activeOverlaps) {
        store.setActiveOverlaps(urlState.activeOverlaps)
      }
      if (urlState.chartOptions) {
        store.updateChartOptions(urlState.chartOptions)
      }
      if (urlState.exNorm !== undefined) {
        store.setExNorm(urlState.exNorm)
      }

      // Mark that store was initialized from URL
      store.setUrlInitialized(true)
    }
  }, [])

  // Fetch and cache spectra metadata
  useSpectraMetadata()

  // Get metadata from store
  const ownerInfo = useOwnerInfo()
  const spectraInfo = useSpectraInfo()

  const [helpOpen, setHelpOpen] = React.useState(false)

  const options = useMemo(() => Object.values(ownerInfo || {}), [ownerInfo])
  const openHelp = useCallback(() => setHelpOpen(true), [])
  const closeHelp = useCallback(() => setHelpOpen(false), [])

  useKeyboardShortcuts()

  // biome-ignore-start format: Keep on one line for ts-expect-error to work
  // @ts-expect-error - WelcomeModal is JSX, will be typed when migrated to TS
  const welcomeModal = <WelcomeModal open={helpOpen} close={closeHelp} isNew={daysSinceLaunch < 120} ownerInfo={ownerInfo} />
  // biome-ignore-end format: End ignore block

  return (
    <>
      <SpectraViewerContainer ownerInfo={ownerInfo} />
      {/* @ts-expect-error - OwnersContainer is JSX, will be typed when migrated to TS */}
      <OwnersContainer ownerInfo={ownerInfo} spectraInfo={spectraInfo} />
      {/* @ts-expect-error - MyAppBar is JSX, will be typed when migrated to TS */}
      <MyAppBar spectraOptions={options} openHelp={openHelp} />
      {welcomeModal}
    </>
  )
}

export default App
