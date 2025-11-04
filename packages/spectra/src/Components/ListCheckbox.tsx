import Checkbox from "@mui/material/Checkbox"
import FormControlLabel from "@mui/material/FormControlLabel"
import MUIListItem from "@mui/material/ListItem"
import type { ChangeEvent } from "react"
import { memo } from "react"

interface ListCheckboxProps {
  /** Callback fired when checkbox state changes */
  onCheckItem: (event: ChangeEvent<HTMLInputElement>) => void
  /** Display name/label for the checkbox */
  name: string
  /** Value associated with the checkbox */
  value: string | number
  /** Whether the checkbox is checked */
  checked: boolean
}

/**
 * A checkbox component rendered as a list item with a form control label.
 * Memoized to prevent unnecessary re-renders.
 */
const ListCheckbox = memo<ListCheckboxProps>(function ListCheckbox({
  onCheckItem,
  name,
  value,
  checked,
}) {
  return (
    <MUIListItem>
      <FormControlLabel
        control={<Checkbox onChange={onCheckItem} checked={checked} value={value} />}
        label={name}
      />
    </MUIListItem>
  )
})

export default ListCheckbox
