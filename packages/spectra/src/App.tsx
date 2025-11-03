import React, { useCallback, useEffect, useMemo } from "react"
import MyAppBar from "./Components/MyAppBar"
import OwnersContainer from "./Components/OwnersContainer"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import useKeyboardShortcuts from "./Components/useKeyboardShortcuts"
import WelcomeModal from "./Components/WelcomeModal"
import { useSpectraMetadata } from "./hooks/useSpectraMetadata"
import { useOwnerInfo, useSpectraInfo } from "./store/metadataStore"
import { useSpectraStore } from "./store/spectraStore"
import { syncURLToStore } from "./store/urlSync"
import "./polyfills"

const daysSinceLaunch = Math.round(
  (new Date(2019, 6, 1).getTime() - Date.now()) / (1000 * 60 * 60 * 24)
)

const App = () => {
  // Sync URL params to store on initial mount
  useEffect(() => {
    syncURLToStore(useSpectraStore.getState())
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

  return (
    <>
      <SpectraViewerContainer ownerInfo={ownerInfo} />
      <OwnersContainer ownerInfo={ownerInfo} spectraInfo={spectraInfo} />
      <MyAppBar spectraOptions={options} openHelp={openHelp} />
      <WelcomeModal
        open={helpOpen}
        close={closeHelp}
        isNew={daysSinceLaunch < 120}
        ownerInfo={ownerInfo}
      />
    </>
  )
}

export default App
