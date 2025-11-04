import CloseIcon from "@mui/icons-material/Close"
import Alert from "@mui/material/Alert"
import Button from "@mui/material/Button"
import IconButton from "@mui/material/IconButton"
import Snackbar from "@mui/material/Snackbar"
import type React from "react"

interface StateConflictToastProps {
  open: boolean
  message: string
  restoreLabel: string
  onRestore: () => void
  onClose: () => void
}

/**
 * Toast notification for state conflict resolution
 * Positioned at top center, 70px from viewport top
 * Dismissible but doesn't auto-hide
 */
export const StateConflictToast: React.FC<StateConflictToastProps> = ({
  open,
  message,
  restoreLabel,
  onRestore,
  onClose,
}) => {
  return (
    <Snackbar
      open={open}
      anchorOrigin={{ vertical: "top", horizontal: "center" }}
      style={{ top: 70 }}
    >
      <Alert
        severity="info"
        action={
          <>
            <Button color="inherit" size="small" onClick={onRestore}>
              {restoreLabel}
            </Button>
            <IconButton size="small" color="inherit" onClick={onClose}>
              <CloseIcon fontSize="small" />
            </IconButton>
          </>
        }
      >
        {message}
      </Alert>
    </Snackbar>
  )
}
