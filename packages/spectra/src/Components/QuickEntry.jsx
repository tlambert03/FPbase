import SearchIcon from "@mui/icons-material/Search"
import IconButton from "@mui/material/IconButton"
import React, { useState } from "react"
import SearchModal from "./SearchModal"

const QuickEntry = ({ options, clearForm }) => {
  const [searchOpen, setSearchOpen] = useState(false)

  const handleClick = () => setSearchOpen(true)

  return (
    <div>
      <IconButton onClick={handleClick}>
        <SearchIcon />
      </IconButton>
      <SearchModal
        options={options}
        open={searchOpen}
        clearForm={clearForm}
        setOpen={setSearchOpen}
      />
    </div>
  )
}

export default QuickEntry
