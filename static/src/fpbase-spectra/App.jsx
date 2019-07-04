import React, { useMemo, useCallback } from "react"
import { SPECTRA_LIST } from "./client/queries"

import { reshapeSpectraInfo } from "./util"
import { SpectraViewer } from "./Components/SpectraViewer"
import OwnersContainer from "./Components/OwnersContainer"
import WelcomeModal from "./Components/WelcomeModal"
import { useCachedQuery } from "./useCachedQuery"
import MyAppBar from "./Components/MyAppBar"
import "./polyfills"

const daysSinceLaunch = Math.round(
  (new Date(2019, 6, 1) - Date.now()) / (1000 * 60 * 60 * 24)
)

const EMPTY = { ownerInfo: {}, spectraInfo: {} }
const App = () => {
  const stash = useCachedQuery(SPECTRA_LIST, "_FPbaseSpectraStash", 5 * 60)
  const { ownerInfo, spectraInfo } = useMemo(() => {
    return stash ? reshapeSpectraInfo(stash.spectra) : EMPTY
  }, [stash])
  window.ownerInfo = ownerInfo
  window.spectraInfo = spectraInfo
  const storageKey = "_hideFPbaseSpectraWelcome"
  const hide = localStorage.getItem(storageKey) === "true"
  const [checked, setChecked] = React.useState(hide)
  const [helpOpen, setHelpOpen] = React.useState(!hide)

  const handleChange = useCallback(e => {
    localStorage.setItem(storageKey, e.target.checked)
    setChecked(e.target.checked)
  }, [])

  const options = useMemo(() => Object.values(ownerInfo || {}), [ownerInfo])
  const openHelp = useCallback(() => setHelpOpen(true), [setHelpOpen])
  const closeHelp = useCallback(() => setHelpOpen(false), [setHelpOpen])
  
  return (
    <>
      <SpectraViewer ownerInfo={ownerInfo} />
      <OwnersContainer
        ownerInfo={ownerInfo}
        spectraInfo={spectraInfo}
      />
      <MyAppBar
        spectraOptions={options}
        openHelp={openHelp}
      />
      <WelcomeModal
        open={helpOpen}
        checked={checked}
        close={closeHelp}
        handleChange={handleChange}
        isNew={daysSinceLaunch < 120}
      />
    </>
  )
}

export default App
