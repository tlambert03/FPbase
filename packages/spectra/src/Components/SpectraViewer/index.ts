import type { FC } from "react"
import type { ChartOptions, OwnerInfo } from "../../types"
import ChartOptionsForm from "./ChartOptionsForm"
import {
  SpectraViewer as SpectraViewerComponent,
  SpectraViewerContainer as SpectraViewerContainerComponent,
} from "./SpectraViewer"

export interface SpectraViewerContainerProps {
  ownerInfo?: Record<string, OwnerInfo>
  provideOptions?: Partial<ChartOptions>
  provideSpectra?: string[]
  provideOverlaps?: string[]
  provideHidden?: string[]
}

// Typed versions
const SpectraViewerContainer: FC<SpectraViewerContainerProps> =
  SpectraViewerContainerComponent as FC<SpectraViewerContainerProps>

// biome-ignore lint/suspicious/noExplicitAny: Legacy component, will be typed properly when migrated
const SpectraViewer: FC<any> = SpectraViewerComponent

export { SpectraViewerContainer, SpectraViewer, ChartOptionsForm }
