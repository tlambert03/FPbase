import React, { useEffect, useRef } from "react"
import Select, { components } from "react-select"
import { categoryIcon } from "./FaIcon"
import WindowedMenuList from "./WindowedMenuList"

// Maximum results to show for short queries (prevents lag on first character)
const MAX_RESULTS_SHORT_QUERY = 500
const SHORT_QUERY_THRESHOLD = 2

const filterOptions = (query, label) => {
  const words = query.toLowerCase().split(" ")
  const opts = label.toLowerCase()
  return words.reduce((acc, cur) => acc && opts.includes(cur), true)
}

const querySorter = (query) => {
  const lowerquery = query.trimRight().toLowerCase()
  return function sortOptions(a, b) {
    const alabel = a.label.replace(/^m/, "").toLowerCase()
    const blabel = b.label.replace(/^m/, "").toLowerCase()
    if (alabel.startsWith(lowerquery)) {
      if (blabel.startsWith(lowerquery)) {
        return alabel < blabel ? -1 : 1
      }
      return -1
    }
    if (blabel.startsWith(lowerquery)) return 1
    return alabel < blabel ? -1 : 1
  }
}

const OptionWithIcon = (props) => {
  const myProps = { ...props }

  myProps.children = (
    <>
      {categoryIcon(myProps.data.category, "#aaa")}
      {myProps.children}
    </>
  )
  return <components.Option {...myProps} />
}

const emptyFilter = () => true

const SortableWindowedSelect = React.memo(function SortableWindowedSelect({
  showIcon = false,
  options,
  components,
  ...otherprops
}) {
  const [dynamicOptions, setOptions] = React.useState(options)
  const [inputValue, setInputValue] = React.useState("")
  const [isLimited, setIsLimited] = React.useState(false)
  const debounceTimerRef = useRef(null)
  const selectRef = useRef()

  useEffect(() => setOptions(options), [options])

  // Escape key blur + cleanup debounce on unmount
  useEffect(() => {
    const handleKeyDown = (e) => e.code === "Escape" && selectRef.current?.blur()
    document.addEventListener("keydown", handleKeyDown)
    return () => {
      document.removeEventListener("keydown", handleKeyDown)
      clearTimeout(debounceTimerRef.current)
    }
  }, [])

  const handleInputChange = React.useCallback(
    (query, { action }) => {
      setInputValue(query)

      if (action === "input-change") {
        clearTimeout(debounceTimerRef.current)
        debounceTimerRef.current = setTimeout(() => {
          if (!query) {
            setOptions(options)
            setIsLimited(false)
            return
          }

          const filtered = options.filter(({ label }) => filterOptions(query, label))
          const sorted = filtered.sort(querySorter(query))
          const shouldLimit =
            query.length < SHORT_QUERY_THRESHOLD && sorted.length > MAX_RESULTS_SHORT_QUERY

          setOptions(shouldLimit ? sorted.slice(0, MAX_RESULTS_SHORT_QUERY) : sorted)
          setIsLimited(shouldLimit)
        }, 150)
      } else if (action === "menu-close") {
        setOptions(options)
      }
    },
    [options]
  )

  const memoizedComponents = React.useMemo(
    () => ({
      ...components,
      ...(showIcon ? { Option: OptionWithIcon } : {}),
      ...{ MenuList: WindowedMenuList },
    }),
    [components, showIcon]
  )

  // Custom message when results are limited
  const noOptionsMessage = React.useCallback(() => {
    if (isLimited) {
      return `Showing first ${MAX_RESULTS_SHORT_QUERY} results. Type more characters to refine search.`
    }
    return "No options"
  }, [isLimited])

  return (
    <Select
      {...otherprops}
      options={dynamicOptions}
      ref={selectRef}
      inputValue={inputValue}
      onInputChange={handleInputChange}
      filterOption={emptyFilter}
      components={memoizedComponents}
      noOptionsMessage={noOptionsMessage}
    />
  )
})

export default SortableWindowedSelect
