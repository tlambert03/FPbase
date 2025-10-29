import React, { useCallback, useMemo } from "react"
import MyAppBar from "./Components/MyAppBar"
import OwnersContainer from "./Components/OwnersContainer"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import WelcomeModal from "./Components/WelcomeModal"
import { useCachedFetch } from "./useCachedQuery"
import { reshapeSpectraInfo } from "./util"
import "./polyfills"
import useKeyboardShortcuts from "./Components/useKeyboardShortcuts"

const daysSinceLaunch = Math.round((new Date(2019, 6, 1) - Date.now()) / (1000 * 60 * 60 * 24))

const EMPTY = { ownerInfo: {}, spectraInfo: {} }
const App = () => {
  const stash = useCachedFetch("/api/proteins/spectraslugs/", "_FPbaseSpectraStash", 60 * 10)
  // const stash = useCachedQuery(SPECTRA_LIST, "_FPbaseSpectraStash", 5 * 60)
  const { ownerInfo, spectraInfo } = useMemo(() => {
    return stash ? reshapeSpectraInfo(stash.spectra) : EMPTY
  }, [stash])
  window.ownerInfo = ownerInfo
  window.spectraInfo = spectraInfo
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
