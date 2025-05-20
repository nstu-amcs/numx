local escape = {
  ["}"] = true,
  ["{"] = true,
  ["&"] = true,
  ["%"] = true,
}

function merge(p, i, j, ...)
  local n = #arg

  for i = 1, n do
    local f = io.open(arg[i], "r")
  end

  local l = 0

  for line in f:lines() do
    l = l + 1

    if l < a then
      goto continue
    end

    if b ~= -1 and l > b then
      goto finish
    end

    local tok = {}

    for t in line:gmatch("[^&]+") do
      tok[#tok + 1] = trim(t)
    end

    for i = 1, #tok, 1 do
      local pfx = ""
      local sfx = ""

      if escape[tok[i]] then
        pfx = "\\"
      end

      if i < #tok then
        sfx = "&"
      end

      tex.print(pfx .. tok[i] .. sfx)
    end

    tex.print("\\\\");

    ::continue::
  end

  ::finish::

  io.close(f)
end
