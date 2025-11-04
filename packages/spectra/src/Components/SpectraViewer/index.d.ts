import type { FC } from "react"
import type { ChartOptions, OwnerInfo } from "../../types"

export interface SpectraViewerContainerProps {
  ownerInfo?: Record<string, OwnerInfo>
  provideOptions?: Partial<ChartOptions>
  provideSpectra?: string[]
  provideOverlaps?: string[]
  provideHidden?: string[]
}

export declare const SpectraViewerContainer: FC<SpectraViewerContainerProps>
// biome-ignore lint/suspicious/noExplicitAny: Legacy component, will be typed properly when migrated
export declare const SpectraViewer: FC<any>
// biome-ignore lint/suspicious/noExplicitAny: Legacy component, will be typed properly when migrated
export declare const ChartOptionsForm: FC<any>
