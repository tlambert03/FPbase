// Type declarations for JavaScript modules that haven't been migrated to TypeScript yet
// These will be removed once all files are migrated to TypeScript

declare module "./Components/OwnersContainer" {
  import type { FC } from "react"
  import type { OwnerInfo } from "../types"

  interface SpectraInfo {
    [spectrumId: string]: {
      subtype: string
      owner: string
      label: string
      category: string
    }
  }

  interface OwnersContainerProps {
    ownerInfo: Record<string, OwnerInfo>
    spectraInfo: SpectraInfo
  }

  const OwnersContainer: FC<OwnersContainerProps>
  export default OwnersContainer
}

declare module "./Components/MyAppBar" {
  import type { FC } from "react"
  import type { OwnerInfo } from "../types"

  interface MyAppBarProps {
    spectraOptions: OwnerInfo[]
    openHelp: () => void
  }

  const MyAppBar: FC<MyAppBarProps>
  export default MyAppBar
}

declare module "./Components/WelcomeModal" {
  import type { FC } from "react"
  import type { OwnerInfo } from "../types"

  interface WelcomeModalProps {
    open: boolean
    close: () => void
    isNew: boolean
    ownerInfo: Record<string, OwnerInfo>
  }

  const WelcomeModal: FC<WelcomeModalProps>
  export default WelcomeModal
}

declare module "./Components/useKeyboardShortcuts" {
  export default function useKeyboardShortcuts(): void
}

// Fallback for other .jsx/.js files
declare module "*.jsx" {
  // biome-ignore lint/suspicious/noExplicitAny: Temporary declaration for unmigrated JSX files
  const Component: React.ComponentType<any>
  export default Component
}

declare module "*.js" {
  // biome-ignore lint/suspicious/noExplicitAny: Temporary declaration for unmigrated JS files
  const module: any
  export = module
}
