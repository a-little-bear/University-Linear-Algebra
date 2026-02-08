/* 内容模式切换：深度连接 ↔ 定义速查 */
;(function () {
  /* ---------- 图标 ---------- */
  var ICON_CONN =
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" width="24" height="24">' +
    '<path fill="currentColor" d="M3.9,12C3.9,10.29 5.29,8.9 7,8.9H11V7H7A5,5 0 0,0 2,12' +
    'A5,5 0 0,0 7,17H11V15.1H7C5.29,15.1 3.9,13.71 3.9,12M8,13H16V11H8V13M17,7H13V8.9H17' +
    'C18.71,8.9 20.1,10.29 20.1,12C20.1,13.71 18.71,15.1 17,15.1H13V17H17A5,5 0 0,0 22,12' +
    'A5,5 0 0,0 17,7Z"/></svg>';

  var ICON_LIST =
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" width="24" height="24">' +
    '<path fill="currentColor" d="M7,5H21V7H7V5M7,13V11H21V13H7M4,4.5A1.5,1.5 0 0,1 5.5,6' +
    'A1.5,1.5 0 0,1 4,7.5A1.5,1.5 0 0,1 2.5,6A1.5,1.5 0 0,1 4,4.5M4,10.5A1.5,1.5 0 0,1 ' +
    '5.5,12A1.5,1.5 0 0,1 4,13.5A1.5,1.5 0 0,1 2.5,12A1.5,1.5 0 0,1 4,10.5M7,19V17H21V19H7' +
    'M4,16.5A1.5,1.5 0 0,1 5.5,18A1.5,1.5 0 0,1 4,19.5A1.5,1.5 0 0,1 2.5,18A1.5,1.5 0 0,1 4,16.5Z"/></svg>';

  var KEY = "la-content-mode";

  /* ---------- 立即应用已保存的模式（减少闪烁） ---------- */
  try {
    if (localStorage.getItem(KEY) === "definitions") {
      document.body.classList.add("definition-sheet");
    }
  } catch (_) {}

  /* ---------- 创建按钮 ---------- */
  function init() {
    var btn = document.createElement("button");
    btn.className = "md-header__button mode-toggle-btn";

    function refresh() {
      var isDef = document.body.classList.contains("definition-sheet");
      btn.innerHTML = isDef ? ICON_LIST : ICON_CONN;
      btn.title = isDef
        ? "当前：定义速查  ·  点击切换为深度连接"
        : "当前：深度连接  ·  点击切换为定义速查";
      btn.setAttribute("aria-label", btn.title);
    }

    refresh();

    btn.addEventListener("click", function () {
      document.body.classList.toggle("definition-sheet");
      var isDef = document.body.classList.contains("definition-sheet");
      try { localStorage.setItem(KEY, isDef ? "definitions" : "connections"); } catch (_) {}
      refresh();
    });

    /* 插入到 palette 切换按钮左侧 */
    var palette = document.querySelector(".md-header__option");
    if (palette && palette.parentNode) {
      palette.parentNode.insertBefore(btn, palette);
    } else {
      /* 降级：放在搜索按钮左侧 */
      var search = document.querySelector(".md-search");
      if (search && search.parentNode) {
        search.parentNode.insertBefore(btn, search);
      }
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
