// Make every content image clickable: open the full-size file in a new tab.
document.addEventListener("DOMContentLoaded", function () {
  document.querySelectorAll(".rst-content img").forEach(function (img) {
    if (img.closest("a")) {
      return; // already linked (e.g., logo or an explicit image link)
    }
    var link = document.createElement("a");
    link.href = img.src;
    link.target = "_blank";
    link.rel = "noopener";
    link.title = "Open full-size image in a new tab";
    img.parentNode.insertBefore(link, img);
    link.appendChild(img);
  });
});
