import React, { useCallback, useEffect, useMemo, useState } from "react"
import MyAppBar from "./Components/MyAppBar"
import OwnersContainer from "./Components/OwnersContainer"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import { StateConflictToast } from "./Components/StateConflictToast"
import useKeyboardShortcuts from "./Components/useKeyboardShortcuts"
import WelcomeModal from "./Components/WelcomeModal"
import { useSpectraMetadata } from "./hooks/useSpectraMetadata"
import { useOwnerInfo, useSpectraInfo } from "./store/metadataStore"
import { useSpectraStore } from "./store/spectraStore"
import type { SpectraURLState } from "./utils/urlParams"
import { canonicalizeURLState, parseURLParams, serializeURLParams } from "./utils/urlParams"
import "./polyfills"

/**
 * Compare two state objects to see if they're meaningfully different
 * Canonicalizes both states (fills in defaults) before serializing for comparison
 */
function statesAreDifferent(stateA: SpectraURLState, stateB: SpectraURLState): boolean {
  const serializedA = serializeURLParams(canonicalizeURLState(stateA))
  const serializedB = serializeURLParams(canonicalizeURLState(stateB))
  return serializedA !== serializedB
}

/**
 * Get the current state from the store as a SpectraURLState object
 */
function getStoreAsURLState(): SpectraURLState {
  const store = useSpectraStore.getState()

  // Only include exNorm if both values are non-null (valid for URL serialization)
  let exNorm: readonly [number, string] | null = null
  if (store.exNorm && store.exNorm[0] !== null && store.exNorm[1] !== null) {
    exNorm = [store.exNorm[0], store.exNorm[1]] as const
  }

  return {
    activeSpectra: store.activeSpectra,
    activeOverlaps: store.activeOverlaps,
    hiddenSpectra: store.hiddenSpectra,
    chartOptions: store.chartOptions,
    exNorm,
    customFilters: store.customFilters,
    customLasers: store.customLasers,
  }
}

/**
 * Detect navigation type using Performance API
 */
function getNavigationType(): string {
  const navEntries = performance.getEntriesByType("navigation")
  if (navEntries.length > 0) {
    const navEntry = navEntries[0] as PerformanceNavigationTiming
    return navEntry.type // 'navigate', 'reload', 'back_forward', 'prerender'
  }
  return "navigate" // fallback
}

const App = () => {
  const [toastState, setToastState] = useState({
    open: false,
    message: "",
    restoreLabel: "",
    restoreState: null as SpectraURLState | null,
    restoreUrl: null as string | null,
  })

  // Handle state initialization with conflict resolution
  useEffect(() => {
    const store = useSpectraStore.getState()
    const urlState = parseURLParams(window.location.search)
    const hasUrlParams = Object.keys(urlState).length > 0

    // Get current store state (loaded from sessionStorage by Zustand persist)
    const sessionState = getStoreAsURLState()
    const hasSessionState =
      (sessionState.activeSpectra?.length ?? 0) > 0 ||
      Object.keys(sessionState.customFilters ?? {}).length > 0 ||
      Object.keys(sessionState.customLasers ?? {}).length > 0

    // No conflict scenarios
    if (!hasUrlParams && !hasSessionState) {
      // Fresh load with no state - do nothing
      return
    }

    if (!hasUrlParams && hasSessionState) {
      // Only session state exists - keep it (already loaded by Zustand)
      return
    }

    if (hasUrlParams && !hasSessionState) {
      // Only URL params exist - apply them
      store.replace(urlState)
      store.setUrlInitialized(true)
      return
    }
    if (!statesAreDifferent(urlState, sessionState)) {
      // States are identical - just mark as URL initialized
      store.setUrlInitialized(true)
      return
    }

    // Conflict exists - resolve based on navigation type
    const navType = getNavigationType()

    if (navType === "reload") {
      // On reload: prioritize session state (keep what user is working on)
      // Session state is already loaded by Zustand, so just clear URL and show toast
      window.history.replaceState({}, "", window.location.pathname)
      setToastState({
        open: true,
        message: "Session no longer matches URL",
        restoreLabel: "Restore from last URL",
        restoreState: urlState,
        restoreUrl: `${window.location.pathname}?${serializeURLParams(urlState)}`,
      })
    } else {
      // On navigate: prioritize URL params (honor explicit navigation)
      store.replace(urlState)
      store.setUrlInitialized(true)

      setToastState({
        open: true,
        message: "Loaded from URL parameters",
        restoreLabel: "Restore last session",
        restoreState: sessionState,
        restoreUrl: null, // Don't change URL when restoring session
      })
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

  // Toast handlers
  const handleRestoreState = useCallback(() => {
    if (toastState.restoreState) {
      const store = useSpectraStore.getState()
      store.replace(toastState.restoreState)

      // If restoring from URL, update the address bar
      if (toastState.restoreUrl) {
        window.history.replaceState({}, "", toastState.restoreUrl)
        store.setUrlInitialized(true)
      }
    }

    setToastState((prev) => ({ ...prev, open: false }))
  }, [toastState.restoreState, toastState.restoreUrl])

  const handleCloseToast = useCallback(() => {
    setToastState((prev) => ({ ...prev, open: false }))
  }, [])

  // biome-ignore-start format: Keep on one line for ts-expect-error to work
  // @ts-expect-error - WelcomeModal is JSX, will be typed when migrated to TS
  const welcomeModal = <WelcomeModal open={helpOpen} close={closeHelp} ownerInfo={ownerInfo} />
  // biome-ignore-end format: End ignore block

  return (
    <>
      <SpectraViewerContainer ownerInfo={ownerInfo} />
      {/* @ts-expect-error - OwnersContainer is JSX, will be typed when migrated to TS */}
      <OwnersContainer ownerInfo={ownerInfo} spectraInfo={spectraInfo} />
      {/* @ts-expect-error - MyAppBar is JSX, will be typed when migrated to TS */}
      <MyAppBar spectraOptions={options} openHelp={openHelp} />
      {welcomeModal}
      <StateConflictToast
        open={toastState.open}
        message={toastState.message}
        restoreLabel={toastState.restoreLabel}
        onRestore={handleRestoreState}
        onClose={handleCloseToast}
      />
    </>
  )
}

export default App
