import React, { memo } from "react"
import Checkbox from "@mui/material/Checkbox"
import MUIListItem from "@mui/material/ListItem"
import FormControlLabel from "@mui/material/FormControlLabel"

const ListCheckbox = memo(function ListCheckbox({
  onCheckItem,
  name,
  value,
  checked,
}) {
  return (
    <MUIListItem>
      <FormControlLabel
        control={
          <Checkbox onChange={onCheckItem} checked={checked} value={value} />
        }
        label={name}
      />
    </MUIListItem>
  )
})

export default ListCheckbox
