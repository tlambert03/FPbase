import React from "react"
import Cookies from "js-cookie"
import { Box, TextField, Button, FormHelperText } from "@mui/material"

function InputForm({ handleSubmit, initialValue }) {

  return (
    <Box component="form" onSubmit={handleSubmit}>
      <input
        type="hidden"
        name="csrfmiddlewaretoken"
        value={Cookies.get("csrftoken")}
      />
      <TextField
        required
        fullWidth
        multiline
        rows={4}
        label="Enter Query"
        name="query"
        defaultValue={initialValue}
        variant="outlined"
        margin="normal"
      />
      <FormHelperText sx={{ mt: -1, mb: 2 }}>
        Single sequence or multiple sequences in FASTA format. Accepts either
        amino acid or nucleotide sequences, but all sequences in a query be of
        the same type.
      </FormHelperText>
      <Button variant="contained" color="primary" type="submit">
        Submit BLAST Query
      </Button>
    </Box>
  )
}

export default InputForm
