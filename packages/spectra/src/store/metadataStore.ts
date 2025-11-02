import { create } from "zustand"
import type { OwnerInfo } from "../types"

/**
 * Metadata store for ownerInfo and spectraInfo
 * Replaces window globals
 */

interface SpectraInfo {
  [spectrumId: string]: {
    subtype: string
    owner: string
    label: string
    category: string
  }
}

interface MetadataState {
  ownerInfo: Record<string, OwnerInfo>
  spectraInfo: SpectraInfo
  setOwnerInfo: (info: Record<string, OwnerInfo>) => void
  setSpectraInfo: (info: SpectraInfo) => void
  setMetadata: (ownerInfo: Record<string, OwnerInfo>, spectraInfo: SpectraInfo) => void
}

export const useMetadataStore = create<MetadataState>((set) => ({
  ownerInfo: {},
  spectraInfo: {},
  setOwnerInfo: (info) => set({ ownerInfo: info }),
  setSpectraInfo: (info) => set({ spectraInfo: info }),
  setMetadata: (ownerInfo, spectraInfo) => set({ ownerInfo, spectraInfo }),
}))

// Selector hooks
export const useOwnerInfo = () => useMetadataStore((state) => state.ownerInfo)
export const useSpectraInfo = () => useMetadataStore((state) => state.spectraInfo)
