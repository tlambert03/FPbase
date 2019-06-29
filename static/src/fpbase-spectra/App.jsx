import React, { useRef, useMemo, useCallback } from "react"
import { SPECTRA_LIST } from "./client/queries"

import { reshapeSpectraInfo } from "./util"
import { SpectraViewer } from "./Components/SpectraViewer"
import OwnersContainer from "./Components/OwnersContainer"
import WelcomeModal from "./Components/WelcomeModal"
import { useCachedQuery } from "./useCachedQuery"
import useSelectors from "./Components/useSelectors"
import MyAppBar from "./Components/MyAppBar"
import memoize from "memoize-one"

const daysSinceLaunch = Math.round(
  (new Date(2019, 6, 1) - Date.now()) / (1000 * 60 * 60 * 24)
)

const memoizedReshapeSpectraInfo = memoize(reshapeSpectraInfo)

const App = () => {
  const stash = useCachedQuery(SPECTRA_LIST, "_FPbaseSpectraStash", 5 * 60)
  const owners = useRef({})
  const spectraInfo = useRef({})
  if (stash) {
    const data = memoizedReshapeSpectraInfo(stash.spectra)
    owners.current = data.owners
    spectraInfo.current = data.spectraInfo
  }

  const {
    selectors,
    changeOwner,
    removeRow,
    clearForm,
    activeSpectra
  } = useSelectors({
    owners: owners.current,
    spectraInfo: spectraInfo.current
  })

  const storageKey = "_hideFPbaseSpectraWelcome"
  const hide = localStorage.getItem(storageKey) === "true"
  const [checked, setChecked] = React.useState(hide)
  const [helpOpen, setHelpOpen] = React.useState(!hide)

  const handleChange = e => {
    localStorage.setItem(storageKey, e.target.checked)
    setChecked(e.target.checked)
  }

  const options = useMemo(() => Object.values(owners.current || {}), [owners])
  const openHelp = useCallback(() => setHelpOpen(true), [setHelpOpen])
  const closeHelp = useCallback(() => setHelpOpen(false), [setHelpOpen])

  return (
    <>
      <SpectraViewer ownerInfo={owners.current} />
      <OwnersContainer
        owners={owners.current}
        selectors={selectors}
        changeOwner={changeOwner}
        removeRow={removeRow}
        clearForm={clearForm}
        activeSpectra={activeSpectra}
      />
      <MyAppBar
        spectraOptions={options}
        clearForm={clearForm}
        openHelp={openHelp}
      />
      <WelcomeModal
        open={helpOpen}
        checked={checked}
        close={closeHelp}
        handleChange={handleChange}
        isNew={daysSinceLaunch < 120}
      />
      {/* <EfficiencyTable
        activeSpectra={activeSpectra}
        spectraInfo={spectraInfo.current}
      /> */}
    </>
  )
}

export default App
